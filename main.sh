#!/bin/bash
# PfCSP ONT pipeline driver. Set the steps to run in config.sh, then: bash main.sh
# Author: Khadim Gueye

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/config.sh"

[[ -f "$CONFIG" ]] || { echo "[ERROR] config not found: $CONFIG" >&2; exit 1; }
# shellcheck source=config.sh
source "$CONFIG"

REF_DIR="${SCRIPT_DIR}/reference"
PROC_DIR="${SCRIPT_DIR}/process"
LOG_DIR="${OUTPUT_DIR}/logs"

YELLOW="\033[33m"; GREEN="\033[32m"; RED="\033[31m"; RESET="\033[0m"
info() { echo -e "${YELLOW}[INFO]${RESET} $*"; }
ok()   { echo -e "${GREEN}[OK]${RESET} $*"; }
die()  { echo -e "${RED}[ERROR]${RESET} $*" >&2; exit 1; }

activate_conda() {
  [[ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]] \
    || die "conda not found at ${CONDA_BASE}"
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}" || die "cannot activate env ${CONDA_ENV}"
}

require_cmd() { command -v "$1" >/dev/null 2>&1 || die "missing command: $1"; }
require_file() { [[ -f "$1" ]] || die "$2 not found: $1"; }
require_dir()  { [[ -d "$1" ]] || die "$2 not found: $1"; }

check_python_modules() {
  python3 - "$@" <<'PY'
import importlib, sys
missing = [m for m in sys.argv[1:] if not importlib.util.find_spec(m)]
if missing:
    raise SystemExit("missing python modules: " + ", ".join(missing))
print("python modules ok: " + ", ".join(sys.argv[1:]))
PY
}

CALLABLE_START=$(( AMPLICON_START + PRIMER_TRIM ))
CALLABLE_END=$(( AMPLICON_END - PRIMER_TRIM ))

case "$MAP_TARGET" in
  genome)   MAP_REF="$GENOME_FASTA" ;;
  amplicon) MAP_REF="${REF_DIR}/pfcsp_amplicon.fasta" ;;
  *) die "MAP_TARGET must be genome or amplicon, got: $MAP_TARGET" ;;
esac

run_step1() {
  info "STEP1 extract CSP reference"
  require_file "$GENOME_FASTA" "genome fasta"
  require_file "$GENOME_GFF" "genome gff"
  mkdir -p "$REF_DIR"
  python3 "${PROC_DIR}/01_extract_reference.py" \
    --genome "$GENOME_FASTA" \
    --gff "$GENOME_GFF" \
    --out-dir "$REF_DIR" \
    --transcript "$CSP_TRANSCRIPT" \
    --gene-id "$CSP_GENE_ID" \
    --amplicon-chrom "$AMPLICON_CHROM" \
    --amplicon-start "$AMPLICON_START" \
    --amplicon-end "$AMPLICON_END" \
    --primer-trim "$PRIMER_TRIM" \
    --th2r "$TH2R_START" "$TH2R_END" \
    --th3r "$TH3R_START" "$TH3R_END"
  ok "reference written to ${REF_DIR}"
}

run_step2() {
  info "STEP2 filter and map reads (target: ${MAP_TARGET})"
  require_dir "$INPUT_DIR" "input dir"
  require_file "$MAP_REF" "mapping reference"
  [[ -f "${MAP_REF}.fai" ]] || samtools faidx "$MAP_REF"
  python3 "${PROC_DIR}/02_map_reads.py" \
    --input-dir "$INPUT_DIR" \
    --out-dir "$OUTPUT_DIR" \
    --reference "$MAP_REF" \
    --chrom "$AMPLICON_CHROM" \
    --start "$AMPLICON_START" \
    --end "$AMPLICON_END" \
    --min-qual "$MIN_READ_QUAL" \
    --min-len "$MIN_READ_LEN" \
    --max-len "$MAX_READ_LEN" \
    --min-mapq "$MIN_MAPQ" \
    --min-overlap "$MIN_READ_OVERLAP" \
    --threads "$THREADS" \
    --jobs "$PARALLEL_JOBS"
  ok "alignments in ${OUTPUT_DIR}/03_alignment"
}

run_step3() {
  info "STEP3 call variants"
  require_file "${REF_DIR}/pfcsp_codon_map.tsv" "codon map (run step 1 first)"
  require_dir "${OUTPUT_DIR}/03_alignment" "alignment dir (run step 2 first)"
  local opt=()
  [[ "${KEEP_INDELS:-0}" -eq 1 ]] && opt+=(--keep-indels)
  if [[ "${USE_CLAIR3:-0}" -eq 1 ]]; then
    require_cmd run_clair3.sh
    opt+=(--use-clair3 --clair3-threads "${CLAIR3_THREADS:-2}")
    [[ -n "${CLAIR3_MODEL:-}" ]] && opt+=(--clair3-model "$CLAIR3_MODEL")
  fi
  python3 "${PROC_DIR}/03_call_variants.py" \
    --out-dir "$OUTPUT_DIR" \
    --reference "$MAP_REF" \
    --codon-map "${REF_DIR}/pfcsp_codon_map.tsv" \
    --chrom "$AMPLICON_CHROM" \
    --callable-start "$CALLABLE_START" \
    --callable-end "$CALLABLE_END" \
    --min-alt-fraction "$MIN_ALT_FRACTION" \
    --min-alt-count "$MIN_ALT_COUNT" \
    --min-base-qual "$MIN_BASE_QUAL" \
    --min-mapq "$MIN_MAPQ" \
    --min-depth "$MIN_SAMPLE_DEPTH" \
    --min-allele-freq "$MIN_ALLELE_FREQ" \
    --min-allele-reads "$MIN_ALLELE_READS" \
    --min-variant-qual "$MIN_VARIANT_QUAL" \
    --strand-filter "$STRAND_FILTER" \
    --ntc-pattern "$NTC_PATTERN" \
    --artefact-background "$ARTEFACT_BACKGROUND" \
    --indel-min-freq "$INDEL_MIN_FREQ" \
    --indel-homopolymer-min "$INDEL_HOMOPOLYMER_MIN" \
    --jobs "$PARALLEL_JOBS" "${opt[@]}"
  ok "variants in ${OUTPUT_DIR}/05_variants"
}

run_step4() {
  info "STEP4 consensus, haplotypes and MOI"
  require_dir "${OUTPUT_DIR}/03_alignment" "alignment dir (run step 2 first)"
  python3 "${PROC_DIR}/04_consensus.py" \
    --out-dir "$OUTPUT_DIR" \
    --chrom "$AMPLICON_CHROM" \
    --callable-start "$CALLABLE_START" \
    --callable-end "$CALLABLE_END" \
    --min-depth "$MIN_SAMPLE_DEPTH" \
    --min-allele-freq "$MIN_ALLELE_FREQ" \
    --min-allele-reads "$MIN_ALLELE_READS" \
    --min-phased-reads "$MIN_PHASED_READS" \
    --min-base-qual "$MIN_BASE_QUAL" \
    --denoise-alpha "$DENOISE_ALPHA" \
    --max-phased-sites "$MAX_PHASED_SITES" \
    --artefact-sites "${OUTPUT_DIR}/05_variants/sites_filtered.tsv" \
    --jobs "$PARALLEL_JOBS"
  ok "consensus in ${OUTPUT_DIR}/06_consensus, haplotypes in ${OUTPUT_DIR}/07_haplotypes"
}

run_step5() {
  info "STEP5 translate and compare with 3D7"
  require_file "${OUTPUT_DIR}/06_consensus/consensus_nt.fasta" "consensus (run step 4 first)"
  local meta_arg=()
  [[ -n "${METADATA:-}" ]] && meta_arg=(--metadata "$METADATA")
  python3 "${PROC_DIR}/05_translate.py" \
    --out-dir "$OUTPUT_DIR" \
    --codon-map "${REF_DIR}/pfcsp_codon_map.tsv" \
    --protein "${REF_DIR}/pfcsp_protein.faa" \
    --callable-start "$CALLABLE_START" \
    --callable-end "$CALLABLE_END" \
    --th2r "$TH2R_START" "$TH2R_END" \
    --th3r "$TH3R_START" "$TH3R_END" \
    "${meta_arg[@]}"
  ok "protein results in ${OUTPUT_DIR}/08_protein"
}

run_step6() {
  info "STEP6 phylogeny against Pf8"
  require_file "${OUTPUT_DIR}/07_haplotypes/haplotypes_nt.fasta" "haplotypes (run step 4 first)"
  require_file "$GENOME_FASTA" "genome fasta"
  require_cmd iqtree
  local pf8=()
  if [[ -n "${PF8_VCF:-}" && -f "${PF8_VCF}" ]]; then
    require_file "$PF8_METADATA" "Pf8 metadata"
    require_file "$PF8_FWS" "Pf8 Fws table"
    pf8=(--pf8-vcf "$PF8_VCF" --pf8-metadata "$PF8_METADATA" --pf8-fws "$PF8_FWS"
         --pf8-focus-countries "$PF8_FOCUS_COUNTRIES"
         --pf8-context-countries "$PF8_CONTEXT_COUNTRIES"
         --pf8-max-focus "$PF8_MAX_FOCUS"
         --pf8-max-context "$PF8_MAX_CONTEXT"
         --pf8-min-fws "$PF8_MIN_FWS"
         --mask-pf8-non-pass "$PF8_MASK_NON_PASS")
  else
    info "no Pf8 VCF at ${PF8_VCF:-unset}, the trees will carry our samples and 3D7 only"
  fi
  python3 "${PROC_DIR}/06_phylogeny.py" \
    --out-dir "$OUTPUT_DIR" \
    --genome "$GENOME_FASTA" \
    --chrom "$AMPLICON_CHROM" \
    --callable-start "$CALLABLE_START" \
    --callable-end "$CALLABLE_END" \
    --sequence-source "$PHYLO_SEQUENCE_SOURCE" \
    --max-ambiguous "$PHYLO_MAX_AMBIGUOUS" \
    --ntc-pattern "$NTC_PATTERN" \
    --positive-pattern "$POSITIVE_PATTERN" \
    --cameroon-prefix "$OUR_CAMEROON_PREFIX" \
    --our-year "$OUR_YEAR" \
    --model "$PHYLO_MODEL" \
    --bootstrap "$PHYLO_BOOTSTRAP" \
    --alrt "$PHYLO_ALRT" \
    --threads "$PHYLO_THREADS" \
    --datasets "$PHYLO_DATASETS" "${pf8[@]}"
  ok "phylogeny in ${OUTPUT_DIR}/09_phylogeny"
}

main() {
  info "config       ${CONFIG}"
  info "input        ${INPUT_DIR}"
  info "output       ${OUTPUT_DIR}"
  info "amplicon     ${AMPLICON_CHROM}:${AMPLICON_START}-${AMPLICON_END}"
  info "callable     ${AMPLICON_CHROM}:${CALLABLE_START}-${CALLABLE_END}"
  info "steps        1:${STEP1} 2:${STEP2} 3:${STEP3} 4:${STEP4} 5:${STEP5} 6:${STEP6}"

  if [[ "${DRY_RUN:-0}" -eq 1 ]]; then
    ok "dry run, nothing executed"
    return 0
  fi

  mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
  activate_conda
  require_cmd python3
  require_cmd minimap2
  require_cmd samtools
  require_cmd bcftools
  require_cmd chopper
  [[ "$STEP3" -eq 1 ]] && require_cmd freebayes
  check_python_modules pysam pandas numpy

  [[ "$STEP1" -eq 1 ]] && run_step1 || info "STEP1 skipped"
  [[ "$STEP2" -eq 1 ]] && run_step2 || info "STEP2 skipped"
  [[ "$STEP3" -eq 1 ]] && run_step3 || info "STEP3 skipped"
  [[ "$STEP4" -eq 1 ]] && run_step4 || info "STEP4 skipped"
  [[ "$STEP5" -eq 1 ]] && run_step5 || info "STEP5 skipped"
  [[ "$STEP6" -eq 1 ]] && run_step6 || info "STEP6 skipped"

  ok "done"
  info "results: ${OUTPUT_DIR}"
  info "figures: run notebook/csp_figures.ipynb"
}

main

# PfCSP ICEMR ONT analysis

Analysis pipeline for the targeted *Plasmodium falciparum* CSP amplicon sequenced on
Oxford Nanopore, covering the ICEMR samples from Nigeria and Cameroon.

The pipeline maps the reads, calls variants without assuming a ploidy, builds a per
sample consensus, resolves the co-circulating clones, and translates the result to
protein in order to score the TH2R and TH3R epitopes against the 3D7 reference.

## Requirements

All tools and versions are pinned in `environment.yml`:

```bash
mamba env create -f environment.yml
```

This creates the environment `pf_ont_pipeline`, which is the value expected by
`CONDA_ENV` in `config.sh`. `main.sh` activates it on its own, so there is no need to
activate it beforehand.

Reference inputs, set in `config.sh`:

- `GENOME_FASTA` Pf3D7 genome, indexed with `samtools faidx`
- `GENOME_GFF` PlasmoDB annotation, used to locate the CSP CDS

## Layout

```
pfcsp_icemr_ont_analysis/
├── environment.yml conda environment, pinned versions
├── config.sh      paths, steps and parameters
├── main.sh        driver, calls the process scripts in order
├── process/       pipeline scripts
├── reference/     CSP reference set, produced by step 1
└── notebook/      figures and report tables
```

## Running

Edit `config.sh`, then:

```bash
bash main.sh
```

`main.sh` takes no arguments. Everything is set in `config.sh`, including which steps
run: `STEP1` to `STEP5`, where 1 runs the step and 0 skips it. Setting `DRY_RUN=1`
prints the resolved configuration and exits without running anything.

Steps 2 to 5 depend on the step before them. Step 1 only needs to be rerun when the
amplicon coordinates or the epitope boundaries change.

## Steps

| Step | Script | Output |
|---|---|---|
| 1 | `01_extract_reference.py` | `reference/` CSP gene, CDS, protein, amplicon, region BED, codon map |
| 2 | `02_map_reads.py` | `02_filtered_reads/`, `03_alignment/`, `01_read_qc/mapping_summary.tsv` |
| 3 | `03_call_variants.py` | `05_variants/` candidate sites, filtered sites, per sample calls |
| 4 | `04_consensus.py` | `06_consensus/`, `07_haplotypes/` consensus, haplotypes, MOI |
| 5 | `05_translate.py` | `08_protein/` protein, substitutions, epitope status |

## Output

```
<OUTPUT_DIR>/
├── 01_read_qc/mapping_summary.tsv        reads, on target fraction, depth per sample
├── 02_filtered_reads/                    reads after chopper
├── 03_alignment/
│   ├── <sample>.bam                      all alignments
│   └── <sample>.ontarget.bam             primary reads overlapping the amplicon
├── 05_variants/
│   ├── per_sample/<sample>.freebayes.vcf.gz
│   ├── sites_filtered.tsv                PASS or ARTEFACT, with amino acid annotation
│   ├── per_sample_calls.tsv              allele frequency and genotype per sample
│   └── allele_depth_matrix.tsv           A C G T counts per sample and site
├── 06_consensus/
│   ├── consensus_nt.fasta                IUPAC consensus per sample
│   └── consensus_positions.tsv           base, depth and ambiguity per position
├── 07_haplotypes/
│   ├── haplotypes.tsv                    one row per resolved clone
│   ├── haplotypes_nt.fasta
│   └── moi_per_sample.tsv                MOI and clone proportions
├── 08_protein/
│   ├── consensus_protein.faa             translated consensus
│   ├── haplotypes_protein.faa            translated clones
│   ├── substitutions_consensus.tsv       substitutions against 3D7, per sample
│   ├── substitutions_haplotypes.tsv      substitutions against 3D7, per clone
│   ├── substitution_prevalence.tsv       prevalence across the cohort
│   ├── epitope_status_consensus.tsv      match, mismatch or incomplete per sample
│   ├── epitope_status_haplotypes.tsv     same, scored per clone
│   └── epitope_alignment_th2r.tsv        residue table for the alignment figure
└── logs/
```

Three files give a first overview: `01_read_qc/mapping_summary.tsv`, then
`07_haplotypes/moi_per_sample.tsv`, then `08_protein/epitope_status_haplotypes.tsv`.

## Parameters

The settings most likely to need changing are grouped at the top of `config.sh`.

| Parameter | Default | Note |
|---|---|---|
| `AMPLICON_START`, `AMPLICON_END` | 221320, 221668 | UNC targeted CSP protocol, 349 bp |
| `PRIMER_TRIM` | 25 | Excluded at each edge, the primer region carries primer sequence |
| `TH2R_START`, `TH2R_END` | 311, 327 | Amino acid positions on the 397 aa protein |
| `TH3R_START`, `TH3R_END` | 352, 367 | Same numbering |
| `MAP_TARGET` | genome | `genome` or `amplicon`, see below |
| `MIN_READ_LEN`, `MAX_READ_LEN` | 250, 600 | Amplicon is 349 bp |
| `MIN_SAMPLE_DEPTH` | 100 | Below this a position is not called |
| `MIN_ALLELE_FREQ` | 0.10 | Minor allele and haplotype floor |
| `MIN_PHASED_READS` | 50 | Complete reads required across the mixed sites |

`MAP_TARGET=genome` maps against the whole Pf3D7 genome and restricts to the amplicon
afterwards. Host and off target reads then land elsewhere instead of being forced onto
a 349 bp reference. `MAP_TARGET=amplicon` maps against the extracted amplicon only,
which is faster but less specific.

## Notes

Amino acid numbering follows PF3D7_0304600.1, 397 residues, so the reference amino acid
printed in each substitution label comes from the 3D7 sequence in `reference/`. The
callable window covers amino acids 292 to 390, which contains both epitopes.

Variants are called with `freebayes --pooled-continuous`. A malaria infection can carry
several clones, so a fixed ploidy model would collapse the allele frequencies that the
clone analysis depends on. freebayes output is used as a candidate list only, and
alleles are recounted from the BAM files at every candidate site in every sample.

A position where two alleles both pass the frequency floor becomes an IUPAC ambiguity
code in the consensus and cannot be scored against 3D7. The haplotype output resolves
this: reads carrying the mixed sites are grouped, denoised with the UNOISE3 rule,
screened for PCR chimeras, and each surviving haplotype is scored on its own.

Control samples are identified from the file name through `NTC_PATTERN` and
`POSITIVE_PATTERN`. NTC depth is used to raise the minor allele floor when the no
template controls are not empty.

## Contact

Khadim Gueye

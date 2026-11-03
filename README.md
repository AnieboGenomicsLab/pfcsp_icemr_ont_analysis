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
│   ├── sites_filtered.tsv                PASS, SINGLETON or ARTEFACT, with annotation
│   ├── per_sample_calls.tsv              allele frequency and genotype per sample
│   ├── allele_depth_matrix.tsv           A C G T counts per sample and site
│   ├── indel_qc.tsv                      indels in the BAM, with homopolymer context
│   └── caller_concordance.tsv            freebayes against Clair3, when enabled
├── 06_consensus/
│   ├── consensus_nt.fasta                IUPAC consensus per sample
│   └── consensus_positions.tsv           base, depth, ambiguity, masked allele
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
| `MAX_PHASED_SITES` | 20 | Above this, the most balanced positions are phased |
| `ARTEFACT_BACKGROUND` | 0.03 | Plate background above which an allele is an artefact |
| `KEEP_INDELS` | 0 | Indels are always reported, 1 also keeps them in the VCF |
| `USE_CLAIR3` | 0 | 1 adds a Clair3 concordance check, changes no result |

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

Each site allele gets one of three statuses. `PASS` is dominant in at least one sample
or a credible minor allele in at least two. `ARTEFACT_recurrent_low_level` is never
dominant anywhere and still sits above `ARTEFACT_BACKGROUND` across a quarter of the
plate, which is the signature of systematic error. `SINGLETON_unconfirmed` is seen in a
single sample: it is kept, because a sample specific allele is the opposite of a
recurrent error, but the label says the evidence rests on one sample.

Only a recurrent artefact is masked before the consensus step picks its mixed
positions, so a systematic error cannot be counted as a clone while a rare clone is not
thrown away. The dominant base of a position is never masked, and
`consensus_positions.tsv` records every allele that was removed.

Phasing uses at most `MAX_PHASED_SITES` positions. Reads span the whole amplicon, so
the limit exists only to keep enough reads clean across every site at once. Above the
limit the most balanced positions are kept rather than the first ones along the genome,
which would drop one epitope in favour of the other. `moi_per_sample.tsv` reports how
many mixed sites were found and how many were used.

Indels are read from the CIGAR strings of every BAM and written to `indel_qc.tsv`
before any filtering. Each one is classified against its homopolymer context, its
frame and its maximum frequency, and anything not explained by homopolymer error
raises a warning. Only substitutions reach the consensus and the protein steps, which
is a safe default for the C terminal window but not for the central repeat region.

`USE_CLAIR3=1` runs Clair3 alongside freebayes and writes `caller_concordance.tsv`.
Clair3 is a deep learning caller trained on ONT R10.4.1, so agreement between the two
shows the calls do not depend on the freebayes error model. Only its allele frequency
is compared, its genotype field assumes a ploidy and is ignored.

Control samples are identified from the file name through `NTC_PATTERN` and
`POSITIVE_PATTERN`. NTC depth is used to raise the minor allele floor when the no
template controls are not empty.

## Contact

Khadim Gueye

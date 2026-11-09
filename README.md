# PfCSP ICEMR ONT analysis

Analysis pipeline for the targeted *Plasmodium falciparum* circumsporozoite protein
(CSP) amplicon sequenced on Oxford Nanopore.

From raw reads it maps and filters, calls variants without assuming a ploidy, builds a
per sample consensus, resolves the co-circulating clones, translates them to protein to
score the Th2R and Th3R epitopes against 3D7, places the sequences against MalariaGEN
Pf8, and draws the publication figure.

The amplicon covers 349 bp of the CSP C terminus, `Pf3D7_03_v3:221320-221668`, which is
amino acids 292 to 390 of `PF3D7_0304600.1` once the primer regions are trimmed. Both
T cell epitopes fall inside it. The central NANP repeat does not.

## Install

Tools and versions are pinned in `environment.yml`:

```bash
mamba env create -f environment.yml
```

This creates the environment `pf_ont_pipeline`, which is what `CONDA_ENV` expects in
`config.sh`. `main.sh` activates it on its own.

## Quick start

Edit the paths at the top of `config.sh`, then:

```bash
bash main.sh
```

Three inputs are required:

| Setting | What it is |
|---|---|
| `INPUT_DIR` | one `<sample>.fastq.gz` per sample |
| `GENOME_FASTA` | Pf3D7 genome, indexed with `samtools faidx` |
| `GENOME_GFF` | PlasmoDB annotation, used to locate the CSP CDS |

The Pf8 comparison in step 6 also needs the MalariaGEN Pf8 release: the CSP region VCF,
`Pf8_samples.txt` and `Pf8_fws.tsv`, set through `PF8_DIR`. Without them step 6 falls
back to your samples and 3D7 alone.

## Steps

Each step is written as `${STEPn:-default}`, so one can be rerun from the command line
without editing the file:

```bash
STEP1=0 STEP2=0 STEP3=0 STEP4=0 STEP5=0 STEP6=1 STEP7=0 bash main.sh
```

| Step | Script | Produces |
|---|---|---|
| 1 | `01_extract_reference.py` | `reference/` CSP gene, CDS, protein, amplicon, codon map |
| 2 | `02_map_reads.py` | `02_filtered_reads/`, `03_alignment/`, mapping summary |
| 3 | `03_call_variants.py` | `05_variants/` candidate sites, filtered sites, per sample calls |
| 4 | `04_consensus.py` | `06_consensus/`, `07_haplotypes/` consensus, clones, MOI |
| 5 | `05_translate.py` | `08_protein/` protein, substitutions, epitope status |
| 6 | `06_phylogeny.py` | `09_phylogeny/` alignments, trees and annotations against Pf8 |
| 7 | `07_figures.py` | `10_figures/` the protein figure against 3D7 |

Steps 2 to 5 depend on the step before them, step 6 needs step 4 and step 7 needs
step 5. Step 1 only has to be rerun when the amplicon coordinates or the epitope
boundaries change. `DRY_RUN=1` prints the resolved configuration and exits.

Step 6 runs IQ-TREE twice and takes about 25 minutes. Run it under `tmux` or `screen`
if the session may be interrupted.

## Output

```
<OUTPUT_DIR>/
├── 01_read_qc/mapping_summary.tsv     reads, on target fraction, depth per sample
├── 03_alignment/<sample>.bam
├── 05_variants/
│   ├── sites_filtered.tsv             PASS, SINGLETON or ARTEFACT, with annotation
│   ├── per_sample_calls.tsv           allele frequency and genotype per sample
│   └── indel_qc.tsv                   indels in the BAM, with homopolymer context
├── 06_consensus/consensus_nt.fasta    IUPAC consensus, one per sample
├── 07_haplotypes/
│   ├── haplotypes_nt.fasta            one sequence per resolved clone
│   └── moi_per_sample.tsv             MOI and clone proportions
├── 08_protein/
│   ├── haplotypes_protein.faa         translated clones
│   ├── substitution_prevalence.tsv    prevalence of each change across the cohort
│   └── epitope_status_haplotypes.tsv  Th2R and Th3R scored per clone
├── 09_phylogeny/<dataset>/
│   ├── <dataset>.rooted.nwk           the tree, rooted on 3D7
│   ├── <dataset>.support_summary.tsv  how much of it the alignment resolves
│   ├── <dataset>.popart_by_group.nex  network input for PopART
│   ├── metadata_sequences.tsv         country, region, year, source per sequence
│   └── itol_*.txt                     annotations to drop on the tree in iTOL
└── 10_figures/
    ├── csp_protein_vs_3d7.pdf         the figure, vector and raster
    └── protein_haplotypes.tsv         every protein haplotype and its changes
```

Three files give a first overview: `01_read_qc/mapping_summary.tsv`, then
`07_haplotypes/moi_per_sample.tsv`, then `08_protein/epitope_status_haplotypes.tsv`.

For the tree, load `<dataset>.rooted.nwk` in iTOL and drop the annotations on it in
this order, each ring label immediately before the ring it names:

```
itol_labelcolors_origin.txt   leaf label colour, adds no ring
itol_ringlabel_origin.txt     itol_colorstrip_origin.txt    ring 1
itol_ringlabel_region.txt     itol_piechart_by_region.txt   ring 2
itol_ringlabel_country.txt    itol_piechart_by_country.txt  ring 3
```

## Configuration

Everything lives in `config.sh`, grouped by stage. The settings most often changed:

| Parameter | Default | Note |
|---|---|---|
| `AMPLICON_START`, `AMPLICON_END` | 221320, 221668 | UNC targeted CSP protocol |
| `PRIMER_TRIM` | 25 | Excluded at each edge of the amplicon |
| `TH2R`, `TH3R` boundaries | 311-327, 352-367 | Amino acid positions on the 397 aa protein |
| `MIN_SAMPLE_DEPTH` | 100 | Below this a position is not called |
| `MIN_ALLELE_FREQ` | 0.10 | Minor allele and haplotype floor |
| `MAX_PHASED_SITES` | 20 | Above this the most balanced positions are phased |
| `MAP_TARGET` | genome | `genome` keeps off target reads off a 349 bp reference |
| `THREADS`, `PARALLEL_JOBS` | 8, 16 | Per job, and jobs at once |

Phylogeny and figure:

| Parameter | Default | Note |
|---|---|---|
| `PHYLO_DATASETS` | nigeria,all | Built in the same run |
| `OUR_CAMEROON_PREFIX` | CMR | Identifier prefix that marks a Cameroonian sample |
| `PF8_FOCUS_COUNTRIES` | Nigeria,Cameroon | Kept deeper than the rest, 80 each |
| `PF8_MAX_CONTEXT` | 20 | Cap for every other country |
| `PF8_MIN_FWS` | 0.95 | Below this a Pf8 sample is a mix of clones |
| `PHYLO_BOOTSTRAP`, `PHYLO_ALRT` | 1000, 1000 | 0 disables either |
| `FIG_MIN_CLONES` | 3 | A haplotype is drawn when this many clones carry it |
| `FIG_DPI` | 600 | Raster resolution of the PNG |

`OUR_CAMEROON_PREFIX` is how the two countries are told apart, from the sample
identifier. Change it before running on another cohort.

## Method notes

The choices that change how a result is read.

**Ploidy is never assumed.** A sample can carry several clones, so a fixed ploidy model
would collapse the allele frequencies the clone analysis depends on. Alleles are
recounted from the BAM at every candidate site in every sample.

**Clones, not consensus.** A polyclonal infection has an IUPAC consensus that no
parasite carries, so the reads are phased instead and the unit of analysis downstream
is one clone. Phasing is n-ary: a position with three alleles above the floor is phased
like any other, and a resolved clone carries a base at every position it was sequenced
at.

**Each site allele gets one of three statuses.** `PASS` is dominant in at least one
sample or a credible minor allele in at least two. `ARTEFACT_recurrent_low_level` is
never dominant and still sits above `ARTEFACT_BACKGROUND` across a quarter of the
plate, which is the signature of systematic error. `SINGLETON_unconfirmed` is seen in
one sample only; it is kept, because a sample specific allele is the opposite of a
recurrent error, but the label says the evidence rests on one sample.

**Indels are reported, not called.** They go to `indel_qc.tsv` with their homopolymer
context and are excluded from the consensus and the protein, which model substitutions
only. In this window the recurrent indel is a 1 bp deletion in a 9 bp poly-T, which is
ONT homopolymer error.

**No aligner runs in step 6.** Every step writes the window in 3D7 coordinates and Pf8
is called against the same reference, so the sequences are already homologous position
by position. Realigning 299 bp of a hypervariable gene would invent gaps the data does
not support.

**The trees are rooted on 3D7 as a reading convention.** 3D7 is a member of the same
species, not an outgroup, and there is no outgroup for *P. falciparum* in this window.

**The window carries little phylogenetic signal.** 22 of its 99 residues vary, and on
the reference cohort that supported 9 percent of internal nodes at `UFBoot >= 95` and
`SH-aLRT >= 80`. Read `<dataset>.support_summary.tsv` before reading the tree, and
prefer the PopART network, which shows the same sequences without committing to an
order of descent the alignment cannot support.

**The figure is drawn at the size it is printed at**, 6.6 by 9.5 inches, the text area
of a portrait page. Insert the PNG at full size in a word processor; scaling it down
shrinks every label by the same factor.

**Scope.** Conclusions from this pipeline concern the CSP C terminal T cell epitopes.
The central NANP repeat, the main B cell target of RTS,S, is outside the amplicon.

## Contact

Khadim Gueye

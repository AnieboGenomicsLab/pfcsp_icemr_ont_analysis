# PfCSP ICEMR ONT analysis

Nanopore amplicon pipeline for the *Plasmodium falciparum* circumsporozoite protein.
From basecalled reads to clone level haplotypes, the protein C terminus against 3D7, and
a phylogeny against MalariaGEN Pf8.

Amplicon `Pf3D7_03_v3:221320-221668`. After trimming 25 bp of primer at each edge the
callable window is 299 bp, amino acids 292 to 390.

## 1. Install

```bash
conda env create -f environment.yml
conda activate pf_ont_pipeline
```

## 2. Configure

Edit `config.sh`. Four paths have to be set, everything else has a working default.

| Setting | What it is |
|---|---|
| `INPUT_DIR` | one `<sample>.fastq.gz` per sample |
| `OUTPUT_DIR` | where results are written |
| `GENOME_FASTA` | Pf3D7 genome, indexed with `samtools faidx` |
| `GENOME_GFF` | PlasmoDB annotation, used to locate the CSP CDS |

For the Pf8 context tree, also set `PF8_VCF`, `PF8_METADATA` and `PF8_FWS`.

## 3. Build the metadata

Once per cohort, by hand. It writes the country and state of each sample.

```bash
python3 process/00_metadata.py \
  --input-dir "$INPUT_DIR" \
  --source /path/to/collection_table.tsv \
  --out "$INPUT_DIR/metadata.tsv"
```

## 4. Run

```bash
DRY_RUN=1 bash main.sh    # print the resolved configuration and stop
bash main.sh
```

To rerun one step, set the others to 0 on the command line.

```bash
STEP1=0 STEP2=0 STEP3=0 STEP4=0 STEP5=0 STEP6=0 STEP7=1 bash main.sh
```

## 5. Draw the figures

```bash
jupyter lab notebook/csp_figures.ipynb
```

Run it top to bottom. Its first cell holds the two paths.

## Steps

| Step | Script | Produces |
|---|---|---|
| 1 | `01_extract_reference.py` | `reference/` gene, CDS, protein, amplicon, codon map |
| 2 | `02_map_reads.py` | `02_filtered_reads/`, `03_alignment/` |
| 3 | `03_call_variants.py` | `05_variants/` candidate sites, filtered sites, per sample calls |
| 4 | `04_consensus.py` | `06_consensus/`, `07_haplotypes/` consensus, clones, MOI |
| 5 | `05_translate.py` | `08_protein/` protein, substitutions, epitope status |
| 6 | `06_phylogeny.py` | `09_phylogeny/` alignments, trees, iTOL annotations |
| 7 | `07_plot_tables.py` | `10_table_for_plot/` one table per figure |

Steps 2 to 5 run in order, step 6 needs step 4, step 7 needs steps 4, 5 and 6. Step 1
only has to be rerun when the amplicon coordinates or the epitope boundaries change.

## Output

```
<OUTPUT_DIR>/
├── 01_read_qc/mapping_summary.tsv     reads, on target fraction, depth per sample
├── 03_alignment/<sample>.bam
├── 05_variants/
│   ├── sites_filtered.tsv             PASS, SINGLETON or ARTEFACT, annotated
│   ├── per_sample_calls.tsv           allele frequency and genotype per sample
│   └── indel_qc.tsv                   indels, with homopolymer context
├── 06_consensus/consensus_nt.fasta    IUPAC consensus, one per sample
├── 07_haplotypes/
│   ├── haplotypes_nt.fasta            one sequence per resolved clone
│   └── moi_per_sample.tsv             MOI and clone proportions
├── 08_protein/
│   ├── haplotypes_protein.faa         translated clones
│   ├── substitution_prevalence.tsv    prevalence of each change
│   └── epitope_status_haplotypes.tsv  Th2R and Th3R scored per clone
├── 09_phylogeny/<dataset>/
│   ├── <dataset>.rooted.nwk           the tree, rooted on 3D7
│   ├── <dataset>.support_summary.tsv  how much of it the alignment resolves
│   └── itol_*.txt                     annotations to load in iTOL
├── 10_table_for_plot/                 step 7, read by the notebook
├── 11_figures/                        the figures
└── 12_tables/                         tables for a manuscript, written by hand
```

Start with `01_read_qc/mapping_summary.tsv`, then `07_haplotypes/moi_per_sample.tsv`,
then `08_protein/epitope_status_haplotypes.tsv`.

## Plot tables

A figure is never redrawn by rerunning the pipeline. Step 7 writes one table per figure
into `10_table_for_plot/`, holding the values and the visual encoding both, colour, size,
marker, order and label. The notebook reads only that folder.

| Figure | Shows | Table |
|---|---|---|
| 1 | samples per country | `fig1_sample_counts.tsv` |
| 2 | MOI per state | `fig2_moi.tsv` |
| 3 | nucleotide and haplotype diversity | `fig3_diversity.tsv` |
| 4 | where the variation sits | `fig4_substitutions.tsv` |
| 7 | vaccine match status | `fig7_vaccine_match.tsv` |
| 8 | FST and PCoA with 3D7 | `fig8_structure.tsv` |
| 9 | protein haplotypes against 3D7 | `fig9_protein_*.tsv`, `fig9_meta.tsv` |
| 10 | the files to load in iTOL | `fig10_tree_annotations.tsv` |

`index.tsv` lists every table and the script that wrote it. To change a colour, edit
`plot_style.tsv` and rerun step 7.

## Phylogeny

Load `<dataset>.rooted.nwk` in iTOL, then drop the annotations on it in this order, each
ring label immediately before the ring it names.

```
itol_labelcolors_origin.txt   leaf label colour, adds no ring
itol_ringlabel_origin.txt     itol_colorstrip_origin.txt    ring 1
itol_ringlabel_our_clones.txt itol_symbol_our_clones.txt    ring 2, circles
itol_ringlabel_region.txt     itol_piechart_by_region.txt   ring 3
itol_ringlabel_country.txt    itol_piechart_by_country.txt  ring 4
```

`itol_bar_our_clones.txt` carries the same clone counts as bars with the number printed
beside each one. Load it instead of the circles when the exact count matters.

## Layout

```
process/                 the pipeline steps, run by main.sh
process_intermediate/    one-off scripts, run by hand, never by main.sh
notebook/                csp_figures.ipynb and csp_plots.py, drawing only
reference/               built by step 1
```

## Settings most often changed

| Parameter | Default | Note |
|---|---|---|
| `AMPLICON_START`, `AMPLICON_END` | 221320, 221668 | UNC targeted CSP protocol |
| `PRIMER_TRIM` | 25 | excluded at each edge |
| `TH2R`, `TH3R` | 311-327, 352-367 | amino acid positions on the 397 aa protein |
| `MIN_SAMPLE_DEPTH` | 100 | below this a position is not called |
| `MIN_ALLELE_FREQ` | 0.10 | minor allele and haplotype floor |
| `MAX_PHASED_SITES` | 20 | above this the most balanced positions are phased |
| `THREADS`, `PARALLEL_JOBS` | 8, 16 | per job, and jobs at once |
| `PHYLO_DATASETS` | nigeria,all | built in the same run |
| `PF8_MIN_FWS` | 0.95 | below this a Pf8 sample is a mix of clones |
| `FIG9_EXCLUDE_PREFIX` | CMR | samples left out of the protein figure |

Each script carries a header explaining what it does and why. Read it before changing a
parameter it depends on.

## Contact

Khadim Gueye

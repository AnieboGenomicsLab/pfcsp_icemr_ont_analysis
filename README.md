# PfCSP ICEMR ONT analysis

Analysis pipeline for the targeted *Plasmodium falciparum* circumsporozoite protein
(CSP) amplicon sequenced on Oxford Nanopore. From raw reads it maps and filters, calls
variants without assuming a ploidy, builds a per sample consensus, resolves the
co-circulating clones, translates them to protein and scores the Th2R and Th3R epitopes
against 3D7, places the sequences against MalariaGEN Pf8, and draws the figure.

This README is a reproduction guide. It states what to install, what to provide, how to
run the pipeline and what comes out.

## 1. Install

Tools and versions are pinned in `environment.yml`:

```bash
mamba env create -f environment.yml
```

This creates the environment `pf_ont_pipeline`, which is what `CONDA_ENV` expects in
`config.sh`. `main.sh` activates it on its own, so the environment does not have to be
active in the shell.

## 2. Provide the inputs

Three are required, set at the top of `config.sh`:

| Setting | What it is |
|---|---|
| `INPUT_DIR` | one `<sample>.fastq.gz` per sample |
| `GENOME_FASTA` | Pf3D7 genome, indexed with `samtools faidx` |
| `GENOME_GFF` | PlasmoDB annotation, used to locate the CSP CDS |

Step 6 also needs the MalariaGEN Pf8 release, set through `PF8_DIR`: the CSP region VCF,
`Pf8_samples.txt` and `Pf8_fws.tsv`. Without them step 6 falls back to your own samples
and 3D7 alone and the rest of the pipeline is unaffected.

`METADATA` is optional. It gives the country and the state of each sample, which the
figure turns into its state columns. Build it once with `00_metadata.py`, which looks
each identifier up in the collection tables of the other studies:

```bash
python3 process/00_metadata.py \
  --input-dir "$INPUT_DIR" \
  --source    /path/to/metadata_pfSmarrt.tsv \
  --source    /path/to/metadata_with_sample.tsv \
  --out       "$INPUT_DIR/metadata.tsv"
```

Any table holding an identifier column and a state column works as a source, and the
first source that holds an identifier wins. The script prints the samples whose state it
had to infer and those it could not resolve. Leave `METADATA` empty to run without it:
the figure then drops its state columns and nothing else changes.

`OUTPUT_DIR` is created if it does not exist. Nothing is written anywhere else.

## 3. Check the configuration

```bash
DRY_RUN=1 bash main.sh
```

This prints the resolved paths, the amplicon window and the steps that would run, then
exits without touching the data. Run it once after editing `config.sh`.

## 4. Run

```bash
bash main.sh
```

Steps 1 and 2 are set to 0 in the shipped `config.sh`, because the reference and the
alignments are built once. On a first run set them to 1.

Step 6 runs IQ-TREE twice and takes about 25 minutes on 16 threads. Run the pipeline
under `tmux` or `screen` if the session may be interrupted.

## 5. Rerun one step

Every step is written as `${STEPn:-default}`, so any of them can be driven from the
command line without editing the file:

```bash
STEP1=0 STEP2=0 STEP3=0 STEP4=0 STEP5=0 STEP6=1 bash main.sh
```

| Step | Script | Produces |
|---|---|---|
| 0 | `00_metadata.py` | `metadata.tsv`, the state of each sample, run once by hand |
| 1 | `01_extract_reference.py` | `reference/` CSP gene, CDS, protein, amplicon, codon map |
| 2 | `02_map_reads.py` | `02_filtered_reads/`, `03_alignment/`, mapping summary |
| 3 | `03_call_variants.py` | `05_variants/` candidate sites, filtered sites, per sample calls |
| 4 | `04_consensus.py` | `06_consensus/`, `07_haplotypes/` consensus, clones, MOI |
| 5 | `05_translate.py` | `08_protein/` protein, substitutions, epitope status |
| 6 | `06_phylogeny.py` | `09_phylogeny/` alignments, trees and annotations against Pf8 |
| 7 | `07_plot_tables.py` | `10_table_for_plot/` one table per figure, plus `11_figures/csp_protein_vs_3d7.pdf` |

Steps 2 to 5 depend on the step before them, step 6 needs step 4, and step 7 needs
steps 4, 5 and 6. Step 1 only has to be rerun when the amplicon coordinates or the
epitope boundaries change.

## 6. Draw the figures

A figure is never redrawn by rerunning the pipeline. Step 7 writes one table per figure
into `10_table_for_plot/`, holding the numbers and the visual encoding both, colour,
size, marker, order and label. The notebook reads only that folder.

```bash
conda activate pf_ont_pipeline
jupyter lab notebook/csp_figures.ipynb
```

Its first cell holds the two paths. Every figure is shown inline and written to
`11_figures/`.

| Figure | What it shows | Table it reads |
|---|---|---|
| 1 | samples per country | `fig1_sample_counts.tsv` |
| 2 | MOI per state | `fig2_moi.tsv` |
| 3 | nucleotide and haplotype diversity | `fig3_diversity.tsv` |
| 4 | where the variation sits | `fig4_substitutions.tsv` |
| 7 | vaccine match status | `fig7_vaccine_match.tsv` |
| 8 | FST and PCoA, with 3D7 on the ordination | `fig8_structure.tsv` |
| 9 | protein haplotypes against 3D7 | `fig9_protein_*.tsv`, `fig9_meta.tsv` |
| 10 | the files to load in iTOL | `fig10_tree_annotations.tsv` |

`index.tsv` lists every table, how many rows it has and which script wrote it. Open it
first when coming back to this a year later.

To change a colour, edit `plot_style.tsv` and rerun step 7. Nothing in the notebook holds
a colour, and nothing in the notebook computes a statistic.

`csp_plots.py` draws and reads tables only. `csp_stats.py`, under `process/`, holds the
population genetics that step 7 runs. `csp_figure9.py` writes the protein matrix tables
and rebuilds that figure from them with `draw_from_tables`. The phylogeny itself is drawn
in iTOL, so the last cell prints the files to load and the order to load them in.

## 7. Output

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
├── 10_table_for_plot/                 written by step 7, read by the notebook
│   ├── index.tsv                      every table, its rows and the script that wrote it
│   ├── plot_style.tsv                 every palette, declared once
│   ├── fig<N>_*.tsv                   one table per figure, values and encoding
│   └── clone_distance_to_3d7.tsv      differences from 3D7 per clone
├── 11_figures/                        written by the notebook and by step 7
│   ├── fig<N>_*.png                   one file per figure
│   └── csp_protein_vs_3d7.pdf         the protein figure, vector and raster
└── 12_tables/                         tables for a manuscript or a report, by hand
```

The three folders are deliberately separate. `10_table_for_plot` feeds the figures and
nothing else. `12_tables` is for the tables a reader sees, which are not the same thing
and should never be read by a plotting script.

Three files give a first overview: `01_read_qc/mapping_summary.tsv`, then
`07_haplotypes/moi_per_sample.tsv`, then `08_protein/epitope_status_haplotypes.tsv`.

For the tree, load `<dataset>.rooted.nwk` in iTOL and drop the annotations on it in
this order, each ring label immediately before the ring it names:

```
itol_labelcolors_origin.txt   leaf label colour, adds no ring
itol_ringlabel_origin.txt     itol_colorstrip_origin.txt    ring 1
itol_ringlabel_our_clones.txt itol_symbol_our_clones.txt    ring 2, circles
itol_ringlabel_region.txt     itol_piechart_by_region.txt   ring 3
itol_ringlabel_country.txt    itol_piechart_by_country.txt  ring 4
```

`itol_bar_our_clones.txt` carries the same clone counts as a bar ring that prints the
number beside each bar. Load it instead of the circles when the exact count matters more
than the shape.

The figure is drawn at the size it is printed at. Insert the PNG at full size in a word
processor; scaling it down shrinks every label by the same factor.

## 8. Settings most often changed

All of them live in `config.sh`, grouped by stage.

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
| `PHYLO_DATASETS` | nigeria,all | Built in the same run |
| `OUR_CAMEROON_PREFIX` | CMR | Identifier prefix that marks a Cameroonian sample |
| `PF8_FOCUS_COUNTRIES` | Nigeria,Cameroon | Kept deeper than the rest, 80 each |
| `PF8_MAX_CONTEXT` | 20 | Cap for every other country |
| `PF8_MIN_FWS` | 0.95 | Below this a Pf8 sample is a mix of clones |
| `PHYLO_BOOTSTRAP`, `PHYLO_ALRT` | 1000, 1000 | 0 disables either |
| `METADATA` | `$INPUT_DIR/metadata.tsv` | Country and state per sample |

`OUR_CAMEROON_PREFIX` is how the two countries are told apart, from the sample
identifier. Change it before running on another cohort.

Each script carries its own header explaining what it does and why it does it that way.
Read that header before changing a parameter it depends on.

## Contact

Khadim Gueye

# Configuration for the PfCSP ONT pipeline.
# Edit paths and parameters here, then run: bash main.sh
# Author: Khadim Gueye

# --- Paths ---

INPUT_DIR="/home/gueyekhadim/work/data/malaria_team/Input_Raw_Data/ICEMR_CRID_CSP"
OUTPUT_DIR="/home/gueyekhadim/work/data/malaria_team/Output_Analysis_Data/ICEMR_CRID_CSP_output"
GENOME_FASTA="/home/gueyekhadim/work/script/malaria_team/pf_haplo_pipeline/reference/Pf3D7_v2.fasta"
GENOME_GFF="/home/gueyekhadim/work/script/malaria_team/pf_haplo_pipeline/reference/PlasmoDB-55_Pfalciparum3D7.gff"

# Optional. Two columns at least: sample_id and country, plus state when available.
# Leave empty to run without metadata.
METADATA=""

CONDA_BASE="/home/gueyekhadim/miniconda3"
CONDA_ENV="pf_ont_pipeline"

# --- Steps: 1 to run, 0 to skip ---

# Set to 1 to print the resolved configuration and exit without running.
DRY_RUN=0

STEP1=0   # extract CSP reference
STEP2=0   # filter and map reads
STEP3=1   # call variants
STEP4=1   # consensus, haplotypes, MOI
STEP5=1   # translate, epitopes, 3D7 comparison

# --- Resources ---

THREADS=8
PARALLEL_JOBS=16

# --- Reference ---

CSP_GENE_ID="PF3D7_0304600"
CSP_TRANSCRIPT="PF3D7_0304600.1"

# Amplicon produced by the UNC targeted CSP protocol (protocols.io hk7zb4zp7).
AMPLICON_CHROM="Pf3D7_03_v3"
AMPLICON_START=221320
AMPLICON_END=221668

# Primer binding width excluded from variant calling at each amplicon edge.
PRIMER_TRIM=25

# Epitope boundaries, amino acid positions on the 397 aa PF3D7_0304600.1 protein.
TH2R_START=311
TH2R_END=327
TH3R_START=352
TH3R_END=367

# --- Read filtering ---

MIN_READ_QUAL=10
MIN_READ_LEN=250
MAX_READ_LEN=600

# --- Mapping ---

# genome   map against the whole Pf3D7 genome, then restrict to the amplicon.
#          Off target and host reads land elsewhere instead of being forced onto CSP.
# amplicon map against the extracted amplicon only. Faster, less specific.
MAP_TARGET="genome"

MIN_MAPQ=20
MIN_READ_OVERLAP=0.8

# --- Variant calling ---

MIN_ALT_FRACTION=0.03
MIN_ALT_COUNT=5
MIN_BASE_QUAL=10

# freebayes QUAL is not a reliable discriminator in pooled-continuous mode, the plate
# level artefact filter does that work. Raise this only to be more conservative.
MIN_VARIANT_QUAL=0

# auto applies the strand balance filter only when both read orientations are present.
# The CSP library is single orientation, so auto skips it. Values: auto, on, off.
STRAND_FILTER="auto"

# --- Indels ---

# Indels are always reported from the BAM in 05_variants/indel_qc.tsv, then excluded
# from the consensus and protein steps, which model substitutions only. In this window
# the single recurrent indel is a 1 bp deletion in a 9 bp poly-T at 221415, found in
# every sample near 10 percent, which is ONT homopolymer error. Set to 1 only if the
# window is extended to the central repeat region, and read the warning it prints.
KEEP_INDELS=0

# Report an indel when it reaches this frequency in at least one sample.
INDEL_MIN_FREQ=0.05

# Homopolymer run length above which a 1 bp indel is treated as basecalling error.
INDEL_HOMOPOLYMER_MIN=5

# --- Optional orthogonal caller ---

# Clair3 is a deep learning caller trained on ONT R10.4.1, run here as a concordance
# check against freebayes. Its genotype field assumes a ploidy and is not used, only
# its allele frequency. This changes no result, it adds caller_concordance.tsv.
USE_CLAIR3=0

# Empty resolves to $CONDA_PREFIX/bin/models/r1041_e82_400bps_sup_v500
CLAIR3_MODEL=""
CLAIR3_THREADS=2

# --- Plate level site filter ---

# An allele is a recurrent ONT artefact when it is never dominant in any sample and
# still sits above this frequency in at least a quarter of the plate. An allele seen in
# one sample only is reported as SINGLETON_unconfirmed instead, and is kept, because a
# sample specific allele is the opposite of a systematic error.
ARTEFACT_BACKGROUND=0.03

# --- Sample and allele thresholds ---

# Minimum depth for a sample to be analysed.
MIN_SAMPLE_DEPTH=100

# Minimum frequency for a minor allele or haplotype to be called.
MIN_ALLELE_FREQ=0.10

# Minimum reads supporting an allele or haplotype.
MIN_ALLELE_READS=10

# Minimum complete reads spanning the phased sites.
MIN_PHASED_READS=50

# UNOISE3 alpha for haplotype denoising.
DENOISE_ALPHA=2.0

# Maximum number of mixed positions used to phase one sample. Reads span the whole
# amplicon, so the limit is only there to keep enough reads clean across every site.
# Above the limit the most balanced positions are kept, since they separate clones best.
MAX_PHASED_SITES=20

# --- Control sample patterns, used for contamination thresholds ---

NTC_PATTERN="NTC"
POSITIVE_PATTERN="3D7|3d7"

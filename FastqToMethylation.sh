#!/usr/bin/env bash
# =============================================================================
# SLURM directives
# =============================================================================
#SBATCH --job-name=Methylation
#SBATCH --output=logs/slurm_%x_%j.out
#SBATCH --error=logs/slurm_%x_%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

# =============================================================================
# FASTQ → Methylation Pipeline (Paired-End)
#
# Reference: Zea mays B73 RefGen_v5
#
# Pipeline:
#   1. Raw FastQC
#   2. Adapter trimming (Trim Galore)
#   3. FastQC of trimmed reads
#   4. Bismark alignment
#   5. Deduplicate aligned reads
#   6. Methylation extraction
#
# Output:
#   • Trimmed FASTQs
#   • Alignment BAM
#   • Deduplicated BAM
#   • CX methylation report
#   • Coverage files
#   • bedGraph
#
# Usage:
#
#   Edit SAMPLE, FASTQ_R1 and FASTQ_R2 below, then submit:
#
#       sbatch FastqToMethylation.sh
#
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# USER CONFIGURATION
# -----------------------------------------------------------------------------

GENOME_FOLDER="/home/jts34805/Genomes/B73V5_Bismark"

OUTDIR="results"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
KEEP_BAMS=TRUE

# ---------------------------------------------------------------------------
# Sample information
# ---------------------------------------------------------------------------

SAMPLE="Leaf_BSseq"

FASTQ_R1="/path/to/sample_R1.fastq.gz"
FASTQ_R2="/path/to/sample_R2.fastq.gz"

# -----------------------------------------------------------------------------
# Load modules
# -----------------------------------------------------------------------------

module purge

module load Bismark
module load Trim_Galore
module load FastQC
module load SAMtools

# -----------------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------------

QC_RAW_DIR="${OUTDIR}/01_fastqc_raw"

TRIM_DIR="${OUTDIR}/02_trimmed"

QC_TRIM_DIR="${OUTDIR}/03_fastqc_trimmed"

ALIGN_DIR="${OUTDIR}/04_alignment"

DEDUP_DIR="${OUTDIR}/05_deduplicated"

METHYL_DIR="${OUTDIR}/06_methylation"

LOG_DIR="${OUTDIR}/logs"

mkdir -p \
    "$QC_RAW_DIR" \
    "$TRIM_DIR" \
    "$QC_TRIM_DIR" \
    "$ALIGN_DIR" \
    "$DEDUP_DIR" \
    "$METHYL_DIR" \
    "$LOG_DIR"

# -----------------------------------------------------------------------------
# Helper: timestamped logging
# -----------------------------------------------------------------------------

log() {

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${SAMPLE}] $*"

}

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------

for var in SAMPLE FASTQ_R1 FASTQ_R2 GENOME_FOLDER
do

    if [[ -z "${!var}" ]]; then

        echo "ERROR: \$${var} is not set."

        exit 1

    fi

done

if [[ ! -f "$FASTQ_R1" ]]; then

    echo "ERROR: FASTQ_R1 not found"

    echo "$FASTQ_R1"

    exit 1

fi

if [[ ! -f "$FASTQ_R2" ]]; then

    echo "ERROR: FASTQ_R2 not found"

    echo "$FASTQ_R2"

    exit 1

fi

if [[ ! -d "$GENOME_FOLDER" ]]; then

    echo "ERROR: Genome folder not found"

    echo "$GENOME_FOLDER"

    exit 1

fi

# -----------------------------------------------------------------------------
# Run information
# -----------------------------------------------------------------------------

log "=============================================================="

log "Starting methylation pipeline"

log "Sample          : ${SAMPLE}"

log "FASTQ R1        : ${FASTQ_R1}"

log "FASTQ R2        : ${FASTQ_R2}"

log "Genome          : ${GENOME_FOLDER}"

log "Threads         : ${THREADS}"

log "Host            : $(hostname)"

log "SLURM Job ID    : ${SLURM_JOB_ID:-not_set}"

log "=============================================================="

# -----------------------------------------------------------------------------
# STEP 1. Raw FastQC
# -----------------------------------------------------------------------------

log "STEP 1: Raw read QC with FastQC"

fastqc \
    -t "$THREADS" \
    -o "$QC_RAW_DIR" \
    "$FASTQ_R1" \
    "$FASTQ_R2" \
    2>"${LOG_DIR}/${SAMPLE}_fastqc_raw.log"

log "Raw FastQC complete."

# -----------------------------------------------------------------------------
# STEP 2. Adapter trimming
# -----------------------------------------------------------------------------

log "STEP 2: Adapter trimming with Trim Galore"

trim_galore \
    --paired \
    --illumina \
    --fastqc \
    --cores 4 \
    --gzip \
    --output_dir "$TRIM_DIR" \
    "$FASTQ_R1" \
    "$FASTQ_R2" \
    2>"${LOG_DIR}/${SAMPLE}_trim_galore.log"

# Trim Galore output filenames

R1_BASE=$(basename "$FASTQ_R1" .fastq.gz)
R2_BASE=$(basename "$FASTQ_R2" .fastq.gz)

TRIM_R1="${TRIM_DIR}/${R1_BASE}_val_1.fq.gz"
TRIM_R2="${TRIM_DIR}/${R2_BASE}_val_2.fq.gz"

if [[ ! -f "$TRIM_R1" || ! -f "$TRIM_R2" ]]; then

    echo "ERROR: Trimmed FASTQ files were not produced."

    exit 1

fi

log "Adapter trimming complete."

log "Trimmed FASTQs:"

log "  ${TRIM_R1}"

log "  ${TRIM_R2}"

# -----------------------------------------------------------------------------
# STEP 3. Alignment
# -----------------------------------------------------------------------------

log "STEP 3: Aligning reads with Bismark"

bismark \
    --bowtie2 \
    --parallel 4 \
    --genome "$GENOME_FOLDER" \
    -1 "$TRIM_R1" \
    -2 "$TRIM_R2" \
    --output_dir "$ALIGN_DIR" \
    2>"${LOG_DIR}/${SAMPLE}_bismark.log"

ALIGN_BAM="${ALIGN_DIR}/$(basename "${TRIM_R1%_val_1.fq.gz}")_bismark_bt2_pe.bam"

if [[ ! -f "$ALIGN_BAM" ]]; then

    echo "ERROR: Bismark alignment failed."

    exit 1

fi

log "Alignment complete."

log "Alignment BAM:"

log "  ${ALIGN_BAM}"

# -----------------------------------------------------------------------------
# STEP 4. Deduplicate aligned reads
# -----------------------------------------------------------------------------

log "STEP 4: Removing PCR duplicates"

deduplicate_bismark \
    --paired \
    --bam \
    --output_dir "$DEDUP_DIR" \
    "$ALIGN_BAM" \
    2>"${LOG_DIR}/${SAMPLE}_deduplicate.log"

DEDUP_BAM="${DEDUP_DIR}/$(basename "${ALIGN_BAM%.bam}").deduplicated.bam"

if [[ ! -f "$DEDUP_BAM" ]]; then

    echo "ERROR: Deduplication failed."

    exit 1

fi

log "Deduplication complete."

log "Deduplicated BAM:"

log "  ${DEDUP_BAM}"
samtools flagstat "$DEDUP_BAM" \
    > "${LOG_DIR}/${SAMPLE}_flagstat.txt"

# -----------------------------------------------------------------------------
# STEP 5. Methylation extraction
# -----------------------------------------------------------------------------

log "STEP 5: Calling methylated cytosines"

bismark_methylation_extractor \
    --paired-end \
    --comprehensive \
    --bedGraph \
    --counts \
    --CX_context \
    --cytosine_report \
    --genome_folder "$GENOME_FOLDER" \
    --gzip \
    --multicore 4 \
    --output "$METHYL_DIR" \
    "$DEDUP_BAM" \
    2>"${LOG_DIR}/${SAMPLE}_methylation_extractor.log"

log "Methylation extraction complete."

# -----------------------------------------------------------------------------
# STEP 6. Locate important output files
# -----------------------------------------------------------------------------

CX_REPORT=$(find "$METHYL_DIR" \
    -name "*.CX_report.txt.gz" \
    | head -n 1)

COV_FILE=$(find "$METHYL_DIR" \
    -name "*.bismark.cov.gz" \
    | head -n 1)

BEDGRAPH=$(find "$METHYL_DIR" \
    -name "*.bedGraph.gz" \
    | head -n 1)

MBIAS=$(find "$METHYL_DIR" \
    -name "*M-bias*.txt" \
    | head -n 1)

# -----------------------------------------------------------------------------
# Optional cleanup
# -----------------------------------------------------------------------------



if [[ "$KEEP_BAMS" != "TRUE" ]]; then

    log "Removing alignment BAM"

    rm -f "$ALIGN_BAM"

fi

# -----------------------------------------------------------------------------
# Finished
# -----------------------------------------------------------------------------

log "=============================================================="

log "Pipeline complete."

log "Output summary"

log ""

log "Raw FastQC:"
log "  ${QC_RAW_DIR}"

log ""

log "Trimmed FASTQs:"
log "  ${TRIM_R1}"
log "  ${TRIM_R2}"

log ""

log "Alignment:"
log "  ${ALIGN_BAM}"

log ""

log "Deduplicated BAM:"
log "  ${DEDUP_BAM}"

log ""

log "CX report:"
log "  ${CX_REPORT}"

log ""

log "Coverage:"
log "  ${COV_FILE}"

log ""

log "bedGraph:"
log "  ${BEDGRAPH}"

if [[ -n "$MBIAS" ]]; then

    log ""

    log "M-bias report:"

    log "  ${MBIAS}"

fi

log "=============================================================="

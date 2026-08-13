#!/usr/bin/env bash
# =============================================================================
# SLURM directives
# =============================================================================
#SBATCH --job-name=BismarkIndex
#SBATCH --output=logs/slurm_%x_%j.out
#SBATCH --error=logs/slurm_%x_%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00

set -euo pipefail

# =============================================================================
# Build Bismark Genome Index
#
# Run this script ONCE for each reference genome.
#
# Input:
#   Directory containing the reference genome FASTA
#
# Output:
#   Bisulfite Bowtie2 index files created inside the genome directory.
#
# =============================================================================

# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

GENOME_FOLDER="/home/jts34805/Genomes/B73V5_Bismark"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# ---------------------------------------------------------------------------
# Load software
# ---------------------------------------------------------------------------

module purge

module load Bismark


# ---------------------------------------------------------------------------
# Logging helper
# ---------------------------------------------------------------------------

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# ---------------------------------------------------------------------------
# Check inputs
# ---------------------------------------------------------------------------

mkdir -p logs

if [[ ! -d "$GENOME_FOLDER" ]]; then
    echo "ERROR: Genome directory not found:"
    echo "  $GENOME_FOLDER"
    exit 1
fi

FASTA_COUNT=$(find "$GENOME_FOLDER" -maxdepth 1 \( -name "*.fa" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fasta.gz" \) | wc -l)

if [[ "$FASTA_COUNT" -eq 0 ]]; then
    echo "ERROR: No FASTA file found in:"
    echo "  $GENOME_FOLDER"
    exit 1
fi

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
log "Reference FASTA(s):"

find "$GENOME_FOLDER" -maxdepth 1 \
    \( -name "*.fa" -o -name "*.fasta" -o -name "*.fa.gz" -o -name "*.fasta.gz" \)

echo
log "==============================================="
log "Building Bismark genome index"
log "Genome directory : $GENOME_FOLDER"
log "Threads          : $THREADS"
log "Host             : $(hostname)"
log "Job ID           : ${SLURM_JOB_ID:-not_set}"
log "==============================================="

# ---------------------------------------------------------------------------
# Build index
# ---------------------------------------------------------------------------

bismark_genome_preparation \
    --bowtie2 \
    --parallel "$THREADS" \
    --verbose \
    "$GENOME_FOLDER"
# ---------------------------------------------------------------------------
# Finished
# ---------------------------------------------------------------------------

log "==============================================="
log "Bismark genome preparation complete."
log "Genome folder:"
log "  $GENOME_FOLDER"
log "==============================================="

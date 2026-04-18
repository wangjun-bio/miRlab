#!/bin/bash
# =============================================================================
# Trim FASTQ Reads to Minimum Length 36nt
# =============================================================================
# This script uses Trimmomatic (single-end mode) to trim reads to a minimum
# length of 36 nucleotides. This step is used to remove short fragments
# (likely degraded RNA or adapter dimers) before microbiome classification.
# =============================================================================

# ----------------------------- Configuration ---------------------------------
# Project root directory (change if needed)
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Input/output directories (relative to PROJECT_ROOT)
INPUT_DIR="${PROJECT_ROOT}/data/raw_fastq"           # Original FASTQ files
OUTPUT_DIR="${PROJECT_ROOT}/data/trimmed_fastq_36nt" # Output trimmed FASTQ files

# Trimmomatic JAR file path (adjust to your system)
# If trimmomatic is in PATH as a command, you can just use "trimmomatic"
TRIMMOMATIC_CMD="trimmomatic"

# Minimum length to keep
MIN_LEN=36

# Number of threads
THREADS=64

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ----------------------------- Helper function ------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ----------------------------- Check prerequisites ---------------------------
if ! command -v "$TRIMMOMATIC_CMD" &> /dev/null; then
    echo "ERROR: Trimmomatic not found. Please install or set TRIMMOMATIC_CMD correctly."
    echo "You can download from: http://www.usadellab.org/cms/?page=trimmomatic"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# ----------------------------- Main processing ------------------------------
# Process all FASTQ files in the input directory
for FILE in "$INPUT_DIR"/*.fastq; do
    # Skip if no files match
    [ -e "$FILE" ] || continue
    
    BASENAME=$(basename "$FILE" .fastq)
    OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}.fastq"
    
    log "Processing $BASENAME ..."
    
    # Run Trimmomatic in single-end mode
    # MINLEN: discard reads shorter than specified length
    $TRIMMOMATIC_CMD SE -threads "$THREADS" \
        "$FILE" "$OUTPUT_FILE" \
        MINLEN:"$MIN_LEN"
    
    if [ $? -eq 0 ]; then
        log "Saved trimmed file to $OUTPUT_FILE"
    else
        log "ERROR: Trimmomatic failed for $BASENAME"
    fi
done

log "All trimming jobs completed."
log "Output directory: $OUTPUT_DIR"
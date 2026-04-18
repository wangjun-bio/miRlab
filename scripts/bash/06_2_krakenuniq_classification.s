#!/bin/bash
# =============================================================================
# KrakenUniq Classification for Unmapped Reads (Microbiome Analysis)
# =============================================================================
# This script:
#   1. Reads trimmed FASTQ files (from 06_1)
#   2. Aligns reads to human genome (T2T-CHM13) using Bowtie2
#   3. Extracts unmapped reads
#   4. Classifies unmapped reads using KrakenUniq against microbial database
#   5. Handles re-sequenced (buce) samples by merging FASTQ files
# =============================================================================

# ----------------------------- Configuration ---------------------------------
# Project root directory (change if needed)
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Input/output directories (relative to PROJECT_ROOT)
TRIMMED_FASTQ_DIR="${PROJECT_ROOT}/data/trimmed_fastq_36nt"   # Input trimmed FASTQ files
OUTPUT_DIR="${PROJECT_ROOT}/results/krakenuniq"               # Output directory
mkdir -p "$OUTPUT_DIR"

# Reference genome (T2T-CHM13) Bowtie2 index
BOWTIE2_INDEX="${PROJECT_ROOT}/data/ref/T2T_CHM13/bowtie2_index/GCF_009914755.1_T2T-CHM13v2.0_genomic"

# KrakenUniq database and executable
KRAKENUNIQ_DB="${PROJECT_ROOT}/data/ref/krakenuniq_db/krakenuniq_standard"
KRAKENUNIQ_CMD="/path/to/krakenuniq"  # Update to your krakenuniq executable path

# Sample list files
BUCE_FILE="${PROJECT_ROOT}/data/sample_lists/buce_202410.txt"

# Threads
THREADS=32

# ----------------------------- Helper function ------------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ----------------------------- Check prerequisites ---------------------------
if [ ! -d "$TRIMMED_FASTQ_DIR" ]; then
    echo "ERROR: Trimmed FASTQ directory not found: $TRIMMED_FASTQ_DIR"
    echo "Please run 06_1_trim_36nt.sh first."
    exit 1
fi

if [ ! -f "$BUCE_FILE" ]; then
    log "Warning: buce file not found: $BUCE_FILE (re-sequenced samples will not be merged)"
fi

# ----------------------------- 1. Read buce mapping file --------------------
declare -A BUCE_MAP
if [ -f "$BUCE_FILE" ]; then
    while IFS=$'\t' read -r orig reseq; do
        BUCE_MAP["$orig"]="$reseq"
    done < "$BUCE_FILE"
fi

# ----------------------------- 2. Process each FASTQ file --------------------
cd "$TRIMMED_FASTQ_DIR" || exit 1

for file in *_cutada_trim_R1.fastq; do
    # Skip if no files match
    [ -e "$file" ] || continue
    
    # Extract base sample ID (format: {sampleID}_cutada_trim_R1.fastq)
    base_name=$(basename "$file" .fastq)
    sample_id=$(echo "$base_name" | awk -F'_' '{print $1"_"$2}')
    log "Processing sample: $sample_id"
    
    # Skip if sample is a re-sequencing file (to avoid double processing)
    if grep -q "$sample_id" <(cut -f2 "$BUCE_FILE" 2>/dev/null); then
        log "Skipping $sample_id (appears as re-sequenced sample, will be merged with original)"
        continue
    fi
    
    # Determine input FASTQ (handle re-sequenced samples)
    input_fastq="$file"
    if [[ -n "${BUCE_MAP[$sample_id]}" ]]; then
        reseq_id="${BUCE_MAP[$sample_id]}"
        reseq_file="${reseq_id}_cutada_trim_R1.fastq"
        if [ -f "$reseq_file" ]; then
            log "Merging original $sample_id with re-sequenced $reseq_id"
            merged_fastq="${sample_id}_merged.fastq"
            cat "$file" "$reseq_file" > "$merged_fastq"
            input_fastq="$merged_fastq"
        else
            log "Warning: re-sequenced file $reseq_file not found, using original only"
        fi
    fi
    
    # ----------------------------- 3. Bowtie2 alignment ------------------------
    # First pass: very sensitive alignment, output unmapped reads to tmp1.fastq
    tmp1_sam="${sample_id}_tmp1.sam"
    tmp1_fastq="${sample_id}_tmp1.fastq"
    log "Bowtie2 pass 1 (very sensitive) for $sample_id"
    bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" \
        -U "$input_fastq" \
        -S "$tmp1_sam" \
        --very-sensitive \
        --un "$tmp1_fastq"
    
    # Second pass: very sensitive local alignment on unmapped reads from pass 1
    tmp2_sam="${sample_id}_tmp2.sam"
    unmapped_fastq="${sample_id}_unmapped.fastq"
    log "Bowtie2 pass 2 (very sensitive local) for $sample_id"
    bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" \
        -U "$tmp1_fastq" \
        -S "$tmp2_sam" \
        --very-sensitive-local \
        --un "$unmapped_fastq"
    
    # Cleanup temporary files
    rm -f "$tmp1_sam" "$tmp2_sam" "$tmp1_fastq"
    if [[ -n "${BUCE_MAP[$sample_id]}" && -f "$merged_fastq" ]]; then
        rm -f "$merged_fastq"
    fi
    
    # ----------------------------- 4. KrakenUniq classification ----------------
    log "KrakenUniq classification for $sample_id"
    report_file="${OUTPUT_DIR}/${sample_id}_reportfile.tsv"
    classification_file="${OUTPUT_DIR}/${sample_id}_readclassification.tsv"
    
    # Preload database (optional, speeds up first run)
    # $KRAKENUNIQ_CMD --db "$KRAKENUNIQ_DB" --preload --threads "$THREADS"
    
    $KRAKENUNIQ_CMD --db "$KRAKENUNIQ_DB" \
        --report-file "$report_file" \
        --threads "$THREADS" \
        "$unmapped_fastq" > "$classification_file"
    
    # Cleanup unmapped FASTQ (can be large)
    rm -f "$unmapped_fastq"
    
    log "Completed $sample_id"
done

# ----------------------------- 5. Cleanup and finalize -----------------------
# Generate file list for reference
ls "$OUTPUT_DIR"/*_reportfile.tsv | sed 's/_reportfile.tsv//' | xargs -n1 basename > "$OUTPUT_DIR/file_list.txt"

log "All KrakenUniq classification jobs completed."
log "Outputs saved to: $OUTPUT_DIR"
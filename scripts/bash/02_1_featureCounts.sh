#!/bin/bash
# =============================================================================
# FeatureCounts for mRNA, lncRNA, snRNA, snoRNA
# =============================================================================
# This script counts reads mapping to transcripts of different RNA types
# using featureCounts (subread package).
# Input:  SAM files from bowtie2 alignment (located in results/alignments/)
# Output: Count tables and reassigned SAM files for each RNA type
# =============================================================================

# ----------------------------- Configuration ---------------------------------
# Project root directory (change if needed)
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Directories (relative to PROJECT_ROOT)
RAW_SAM_DIR="${PROJECT_ROOT}/results/alignments"      # input SAM files
OUTPUT_DIR="${PROJECT_ROOT}/results/featureCounts"    # output directory
REF_DIR="${PROJECT_ROOT}/data/ref/annotation"         # GTF annotation files

# Sample list file (one sample ID per line)
SAMPLE_LIST="${PROJECT_ROOT}/data/sample_lists/split_barcode_reads_list.txt"

# GTF files (relative to REF_DIR)
GTF_MRNA="${REF_DIR}/Coding_gene_annotation.gtf"
GTF_LNCRNA="${REF_DIR}/lncRNA_noMIR_noSNHG_noRNY_annotation.gtf"
GTF_SNRNA="${REF_DIR}/snRNA_annotation.gtf"
GTF_SNORNA="${REF_DIR}/snoRNA_annotation.gtf"

# Number of threads for featureCounts
THREADS=60

# ----------------------------- Check prerequisites ---------------------------
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found: $SAMPLE_LIST"
    exit 1
fi

for gtf in "$GTF_MRNA" "$GTF_LNCRNA" "$GTF_SNRNA" "$GTF_SNORNA"; do
    if [ ! -f "$gtf" ]; then
        echo "ERROR: GTF file not found: $gtf"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

# ----------------------------- Main loop -------------------------------------
while read -r sample_id; do
    # Skip empty lines
    [ -z "$sample_id" ] && continue

    echo "Processing sample: $sample_id"

    # Input SAM file (from bowtie2 alignment)
    input_sam="${RAW_SAM_DIR}/${sample_id}mapping2transcript.sam"
    if [ ! -f "$input_sam" ]; then
        echo "WARNING: Input SAM not found: $input_sam - skipping"
        continue
    fi

    # ----- mRNA -----
    featureCounts -T "$THREADS" -O -M -R SAM -s 1 \
        -t transcript -g gene_name \
        -a "$GTF_MRNA" \
        -o "${OUTPUT_DIR}/${sample_id}mRNA_S1_transcript_cal.txt" \
        "$input_sam" \
        2> "${OUTPUT_DIR}/${sample_id}mRNA_S1_transcript_cal.log"
    # Rename the reassigned SAM file (featureCounts adds .featureCounts.sam suffix)
    if [ -f "${input_sam}.featureCounts.sam" ]; then
        mv "${input_sam}.featureCounts.sam" "${OUTPUT_DIR}/${sample_id}mRNA_S1_transcript.sam"
    fi

    # ----- lncRNA -----
    featureCounts -T "$THREADS" -O -M -R SAM -s 1 \
        -t transcript -g gene_name \
        -a "$GTF_LNCRNA" \
        -o "${OUTPUT_DIR}/${sample_id}lncRNA_S1_transcript_cal.txt" \
        "$input_sam" \
        2> "${OUTPUT_DIR}/${sample_id}lncRNA_S1_transcript_cal.log"
    if [ -f "${input_sam}.featureCounts.sam" ]; then
        mv "${input_sam}.featureCounts.sam" "${OUTPUT_DIR}/${sample_id}lncRNA_S1_transcript.sam"
    fi

    # ----- snRNA -----
    featureCounts -T "$THREADS" -O -M -R SAM -s 1 \
        -t transcript -g gene_name \
        -a "$GTF_SNRNA" \
        -o "${OUTPUT_DIR}/${sample_id}snRNA_S1_transcript_cal.txt" \
        "$input_sam" \
        2> "${OUTPUT_DIR}/${sample_id}snRNA_S1_transcript_cal.log"
    if [ -f "${input_sam}.featureCounts.sam" ]; then
        mv "${input_sam}.featureCounts.sam" "${OUTPUT_DIR}/${sample_id}snRNA_S1_transcript.sam"
    fi

    # ----- snoRNA -----
    featureCounts -T "$THREADS" -O -M -R SAM -s 1 \
        -t transcript -g gene_name \
        -a "$GTF_SNORNA" \
        -o "${OUTPUT_DIR}/${sample_id}snoRNA_S1_transcript_cal.txt" \
        "$input_sam" \
        2> "${OUTPUT_DIR}/${sample_id}snoRNA_S1_transcript_cal.log"
    if [ -f "${input_sam}.featureCounts.sam" ]; then
        mv "${input_sam}.featureCounts.sam" "${OUTPUT_DIR}/${sample_id}snoRNA_S1_transcript.sam"
    fi

done < "$SAMPLE_LIST"

echo "All featureCounts jobs completed. Results in: $OUTPUT_DIR"
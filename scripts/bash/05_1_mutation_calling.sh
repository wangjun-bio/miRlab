#!/bin/bash
# =============================================================================
# Mutation Calling Pipeline (Bowtie2 -> SAMtools -> FreeBayes -> ANNOVAR)
# =============================================================================
# This script:
#   1. Processes paired-end FASTQ files (after adapter trimming)
#   2. Handles re-sequenced samples (buce) by merging reads
#   3. Aligns reads to human genome (hg38) using Bowtie2
#   4. Filters for uniquely mapped reads (no secondary alignments, no XS tag)
#   5. Converts SAM to BAM, sorts, and indexes
#   6. Calls variants using FreeBayes
#   7. Filters variants by depth, allele frequency, and quality
#   8. Annotates variants using ANNOVAR
# =============================================================================

# ----------------------------- Configuration ---------------------------------
# Project root directory (change if needed)
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Input/output directories (relative to PROJECT_ROOT)
RAW_FASTQ_DIR="${PROJECT_ROOT}/data/raw_fastq"           # input FASTQ files
BOWTIE2_INDEX="${PROJECT_ROOT}/data/ref/hg38/bowtie2_index/hg38"  # Bowtie2 index prefix
REF_FASTA="${PROJECT_ROOT}/data/ref/hg38/hg38.fa"        # Reference genome FASTA
BAM_OUT_DIR="${PROJECT_ROOT}/results/mutation/bam"       # BAM files
VCF_OUT_DIR="${PROJECT_ROOT}/results/mutation/vcf"       # VCF files
ANNOVAR_OUT_DIR="${PROJECT_ROOT}/results/mutation/annovar" # ANNOVAR outputs

# Sample list files (relative to PROJECT_ROOT)
BUCE_FILE="${PROJECT_ROOT}/data/sample_lists/buce_202410.txt"  # re-sequenced sample mapping
SAMPLE_LIST="${PROJECT_ROOT}/data/sample_lists/sample_list.txt"  # will be generated

# ANNOVAR database and scripts (adjust paths as needed)
ANNOVAR_DB="/path/to/annovar/humandb_hg38"   # Update to your ANNOVAR database path
ANNOVAR_SCRIPT="/path/to/annovar/table_annovar.pl"  # Update to your ANNOVAR script path
CONVERT2ANNOVAR="/path/to/annovar/convert2annovar.pl"

# Threads
THREADS=60

# Create output directories
mkdir -p "$BAM_OUT_DIR" "$VCF_OUT_DIR" "$ANNOVAR_OUT_DIR"

# ----------------------------- Helper functions -----------------------------
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# ----------------------------- 1. Prepare sample list -----------------------
# Generate list of sample IDs from FASTQ files
# Expected FASTQ naming: {sampleID}_cutada_trim_R1.fastq
cd "$RAW_FASTQ_DIR" || exit 1
ls *_cutada_trim_R1.fastq | sed 's/_cutada_trim_R1.fastq//' > "$PROJECT_ROOT/temp_sample_list.txt"

# ----------------------------- 2. Read buce mapping file --------------------
# buce_202410.txt format: original_id<TAB>re_seq_id
declare -A BUCE_MAP
while IFS=$'\t' read -r orig reseq; do
    BUCE_MAP["$orig"]="$reseq"
done < "$BUCE_FILE"

# ----------------------------- 3. Process each sample -----------------------
while read -r sample_id; do
    log "Processing sample: $sample_id"
    
    # Skip if sample is a re-sequencing file (to avoid double processing)
    # Only process original IDs (those not appearing as reseq in buce file)
    if grep -q "$sample_id" <(cut -f2 "$BUCE_FILE"); then
        log "Skipping $sample_id (appears as re-sequenced sample, will be merged with original)"
        continue
    fi
    
    # Determine input FASTQ
    fastq_file="${RAW_FASTQ_DIR}/${sample_id}_cutada_trim_R1.fastq"
    if [[ -n "${BUCE_MAP[$sample_id]}" ]]; then
        reseq_id="${BUCE_MAP[$sample_id]}"
        reseq_file="${RAW_FASTQ_DIR}/${reseq_id}_cutada_trim_R1.fastq"
        if [[ -f "$reseq_file" ]]; then
            log "Merging original $sample_id with re-sequenced $reseq_id"
            merged_fastq="${RAW_FASTQ_DIR}/${sample_id}_merged.fastq"
            cat "$fastq_file" "$reseq_file" > "$merged_fastq"
            input_fastq="$merged_fastq"
        else
            log "Warning: re-sequenced file $reseq_file not found, using original only"
            input_fastq="$fastq_file"
        fi
    else
        input_fastq="$fastq_file"
    fi
    
    # Bowtie2 alignment
    sam_file="${sample_id}_bowtie2.sam"
    log "Bowtie2 alignment for $sample_id"
    bowtie2 -p "$THREADS" -x "$BOWTIE2_INDEX" -U "$input_fastq" -S "$sam_file"
    
    # Filter for uniquely mapped reads (no secondary, no XS tag)
    unique_sam="${sample_id}_bowtie2_unique.sam"
    log "Filtering unique reads for $sample_id"
    samtools view -h "$sam_file" | awk '$0 ~ /^@/ || ($2 !~ /4/ && $0 !~ /XS:i:/)' > "$unique_sam"
    
    # Convert to BAM, sort, index
    unique_bam="${sample_id}_bowtie2_unique.bam"
    sorted_bam="${sample_id}_bowtie2_unique_sorted.bam"
    log "Converting to BAM and sorting for $sample_id"
    samtools view -Sb "$unique_sam" > "$unique_bam"
    samtools sort -o "$sorted_bam" "$unique_bam"
    samtools index "$sorted_bam"
    
    # Move BAM files to output directory
    mv "$sorted_bam" "$BAM_OUT_DIR/"
    mv "${sorted_bam}.bai" "$BAM_OUT_DIR/"
    
    # Cleanup intermediate files
    rm -f "$sam_file" "$unique_sam" "$unique_bam"
    if [[ -n "${BUCE_MAP[$sample_id]}" && -f "$merged_fastq" ]]; then
        rm -f "$merged_fastq"
    fi
    
    # Record sample ID for later steps
    echo "$sample_id" >> "$PROJECT_ROOT/bam_list.txt"
    
done < "$PROJECT_ROOT/temp_sample_list.txt"

rm -f "$PROJECT_ROOT/temp_sample_list.txt"

# ----------------------------- 4. Variant calling with FreeBayes ------------
cd "$BAM_OUT_DIR" || exit 1
BAM_LIST="$PROJECT_ROOT/bam_list.txt"

# Run FreeBayes in parallel (using GNU parallel)
cat "$BAM_LIST" | parallel -j "$THREADS" '
    sample={}
    log "FreeBayes for $sample"
    freebayes -f "'"$REF_FASTA"'" -p 1 --min-base-quality 20 \
        ./${sample}_bowtie2_unique_sorted.bam > "'"$VCF_OUT_DIR"'"/${sample}.vcf
'

# ----------------------------- 5. Filter VCF files --------------------------
cd "$VCF_OUT_DIR" || exit 1
cat "$BAM_LIST" | parallel -j "$THREADS" '
    sample={}
    log "Filtering VCF for $sample"
    bcftools view -i "FORMAT/DP>=30 && FORMAT/AO/(FORMAT/AO+FORMAT/RO)>=0.4 && QUAL>=20" \
        ${sample}.vcf -o ${sample}_filtered.vcf
'

# ----------------------------- 6. ANNOVAR annotation ------------------------
cd "$VCF_OUT_DIR" || exit 1
cat "$BAM_LIST" | parallel -j "$THREADS" '
    sample={}
    log "ANNOVAR annotation for $sample"
    perl "'"$CONVERT2ANNOVAR"'" -format vcf4old ${sample}_filtered.vcf \
        -outfile "'"$ANNOVAR_OUT_DIR"'"/${sample}_filtered.avinput
    perl "'"$ANNOVAR_SCRIPT"'" "'"$ANNOVAR_OUT_DIR"'"/${sample}_filtered.avinput \
        "'"$ANNOVAR_DB"'" -buildver hg38 -out "'"$ANNOVAR_OUT_DIR"'"/${sample} \
        -protocol refGeneWithVer,cytoBand,gnomad211_exome,avsnp151,dbnsfp47a \
        -operation gx,r,f,f,f -xref /path/to/gene_xref.txt -remove -nastring . -csvout -polish
'

# ----------------------------- 7. Move final outputs ------------------------
# BAM files already in BAM_OUT_DIR
# VCF and ANNOVAR outputs remain in their respective directories
mv "$BAM_LIST" "$ANNOVAR_OUT_DIR/"

log "Mutation calling pipeline completed."
log "Outputs:"
log "  BAM files: $BAM_OUT_DIR"
log "  VCF files: $VCF_OUT_DIR"
log "  ANNOVAR outputs: $ANNOVAR_OUT_DIR"
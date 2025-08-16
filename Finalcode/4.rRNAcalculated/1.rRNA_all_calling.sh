# This script is used to mapping all reads to silva rRNA database
# modification : 25-08-11

###1.  generate the star index for silva rRNA database
# SSU
STAR --runMode genomeGenerate \
     --genomeDir /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_NR99_STAR \
     --genomeFastaFiles /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
     --runThreadN 60 \
     --sjdbOverhang 35 \
     --limitGenomeGenerateRAM 68604848736
     

# LSU
STAR --runMode genomeGenerate \
    --genomeDir /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_NR99_STAR_LSU \
    --genomeFastaFiles /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_138.2_LSURef_NR99_tax_silva.fasta \
    --runThreadN 60 \
    --limitGenomeGenerateRAM 68604848736 \
    --sjdbOverhang 35

### 2. combine the all sample index 
cd /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region
mkdir -p rRNA_calculate
rm ./all_samples.txt
find /mnt/data3/yiyonghao/NC_paper_rawdata/upload/star_mapping/unmapped_reads/*.fastq -type f >> ./all_samples.txt
find /mnt/data3/yiyonghao/NC_paper_rawdata/PN/star_mapping/unmapped_reads/*.fastq -type f >> ./all_samples.txt
find /mnt/data3/yiyonghao/NC_paper_rawdata/Brain/unmapped_reads/*.fastq -type f >> ./all_samples.txt


# then use the index to map the reads to silva rRNA database
#!/bin/bash
# set -euo pipefail
IFS=$'\n\t'

echo "Start processing: $(date)"

cd /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_calculate

# 1) initialization folds and file
mkdir -p SILVA_NR99_STAR_SSU LSU_reads SILVA_NR99_STAR_LSU
: > mapping_summary_SSU.txt
: > mapping_summary_LSU.txt

# 2) Traverse the samples 
while IFS= read -r sample; do
  [[ -z "$sample" ]] && continue
  filename=$(basename "$sample")
  echo ">>> Sample: $filename"

  # --- SSU mapping ---
  STAR \
    --genomeDir /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_NR99_STAR \
    --readFilesIn "$sample" \
    --runThreadN 60 \
    --outFileNamePrefix "SILVA_NR99_STAR_SSU/${filename}_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --seedSearchStartLmax 18 \
    --seedPerReadNmax 1000000 \
    --seedPerWindowNmax 50 \
    --outFilterMismatchNmax 3 \
    --outFilterMismatchNoverLmax 0.12 \
    --outFilterMatchNminOverLread 0.20 \
    --outFilterScoreMinOverLread 0.20 \
    --alignIntronMax 1 \
    --alignMatesGapMax 0 \
    --chimOutType None \
    --outFilterMultimapNmax 100

    # 3) Extract the unmapped reads and trans to FASTQ
    samtools fastq -f 4 -F 0x100 -F 0x800 "$ssu_bam" > "LSU_reads/${filename}"

    length=$(wc -l LSU_reads/${filename})
    echo "The reads num of ${filename} not mapped to SSU is ${length}"
    --- LSU mapping ---
    STAR \
      --genomeDir /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_NR99_STAR_LSU \
      --readFilesIn "LSU_reads/${filename}" \
      --runThreadN 60 \
      --outFileNamePrefix "SILVA_NR99_STAR_LSU/${filename}_" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMunmapped Within \
      --seedSearchStartLmax 18 \
      --seedPerReadNmax 1000000 \
      --seedPerWindowNmax 50 \
      --outFilterMismatchNmax 3 \
      --outFilterMismatchNoverLmax 0.12 \
      --outFilterMatchNminOverLread 0.20 \
      --outFilterScoreMinOverLread 0.20 \
      --alignIntronMax 1 \
      --alignMatesGapMax 0 \
      --chimOutType None \
      --outFilterMultimapNmax 100


done < "../all_samples.txt"

echo "Done: $(date)"

#!/usr/bin/env bash
set -euo pipefail

JOBS=60

# output files 
SSU_SUMMARY="mapping_summary_SSU.txt"
LSU_SUMMARY="mapping_summary_LSU.txt"
LOCK_SSU=".lock_ssu"
LOCK_LSU=".lock_lsu"

# initialization
printf "sample\tmapped\tunmapped\n" > "$SSU_SUMMARY"
printf "sample\tmapped\tunmapped\n" > "$LSU_SUMMARY"

process_one_sample() {
  local sample="$1"
  [[ -z "$sample" ]] && return 0
  local filename
  filename=$(basename "$sample")
  echo ">>> Sample: $filename"

  local ssu_bam="SILVA_NR99_STAR_SSU/${filename}_Aligned.sortedByCoord.out.bam"
  if [[ -f "$ssu_bam" ]]; then
    local mapped unmapped
    mapped=$(samtools view -c -F 4 -F 0x100 -F 0x800 "$ssu_bam")
    unmapped=$(samtools view -c -f 4 -F 0x100 -F 0x800 "$ssu_bam")
    { flock -x 200; printf "%s\t%s\t%s\n" "$filename" "$mapped" "$unmapped" >> "$SSU_SUMMARY"; } 200>"$LOCK_SSU"

    # --- LSU ---
    local lsu_bam="SILVA_NR99_STAR_LSU/${filename}_Aligned.sortedByCoord.out.bam"
    if [[ -f "$lsu_bam" ]]; then
      local mapped2 unmapped2
      mapped2=$(samtools view -c -F 4 -F 0x100 -F 0x800 "$lsu_bam")
      unmapped2=$(samtools view -c -f 4 -F 0x100 -F 0x800 "$lsu_bam")
      { flock -x 201; printf "%s\t%s\t%s\n" "$filename" "$mapped2" "$unmapped2" >> "$LSU_SUMMARY"; } 201>"$LOCK_LSU"
    else
      echo "Warning: LSU BAM not found for $filename" >&2
      { flock -x 201; printf "%s\tNA\tNA\n" "$filename" >> "$LSU_SUMMARY"; } 201>"$LOCK_LSU"
    fi
  else
    echo "Warning: SSU BAM not found for $filename" >&2
    { flock -x 200; printf "%s\tNA\tNA\n" "$filename" >> "$SSU_SUMMARY"; } 200>"$LOCK_SSU"
  fi
}

export -f process_one_sample
export SSU_SUMMARY LSU_SUMMARY LOCK_SSU LOCK_LSU

# parallel processed
parallel --jobs "$JOBS" --linebuffer process_one_sample :::: ../all_samples.txt

# This script is used to mapping all reads to silva rRNA database
# modification : 25-08-11

#!/usr/bin/env bash
IFS=$'\n\t'

# ========= 0) initialized =========
ROOT="/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region"
SSU_IDX="${ROOT}/SILVA_NR99_STAR"
LSU_IDX="${ROOT}/SILVA_NR99_STAR_LSU"
LIST_ALL="${ROOT}/all_samples_C.txt"
WORK="${ROOT}/rRNA_calculate"
SUBSET_DIR="${WORK}/20k_subset"
LIST_20K="${ROOT}/all_samples_20k.txt"

# ========= 1) Generate index =========
# STAR --runMode genomeGenerate \
#   --genomeDir "${SSU_IDX}" \
#   --genomeFastaFiles /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
#   --runThreadN 60 \
#   --limitGenomeGenerateRAM 68604848736
#
# STAR --runMode genomeGenerate \
#   --genomeDir "${LSU_IDX}" \
#   --genomeFastaFiles /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/SILVA_138.2_LSURef_NR99_tax_silva.fasta \
#   --runThreadN 60 \
#   --limitGenomeGenerateRAM 68604848736

# ========= 2) Find all samples =========
cd "${ROOT}"
: > "${LIST_ALL}"
find /mnt/data3/yiyonghao/NC_paper_rawdata/upload/classification_reads  -type f -name '*_classRead.fastq' -printf '%p\n' >> "${LIST_ALL}"
find /mnt/data3/yiyonghao/NC_paper_rawdata/PN/classification_reads    -type f -name '*_classRead.fastq' -printf '%p\n' >> "${LIST_ALL}"
find /mnt/data3/yiyonghao/NC_paper_rawdata/Brain/classification_reads            -type f -name '*_classRead.fastq' -printf '%p\n' >> "${LIST_ALL}"

# ========= 3) Random select 20k reads from each sample =========
mkdir -p "${SUBSET_DIR}"

while IFS= read -r sample; do
  [[ -z "${sample}" ]] && continue
  name="$(basename "${sample}" _classRead.fastq)"
  out="${SUBSET_DIR}/${name}_classRead20k.fastq"
  echo "[SUBSET] ${sample} -> ${out}"
  seqtk sample -s42 "${sample}" 20000 > "${out}"
done < "${LIST_ALL}"

# Find subset samples 
: > "${LIST_20K}"
find "${SUBSET_DIR}" -type f -name '*_classRead20k.fastq' -printf '%p\n' > "${LIST_20K}"

# ========= 4) SSU mapping & LSU mapping =========
cd "${SUBSET_DIR}"
mkdir -p SILVA_NR99_STAR_SSU SILVA_NR99_STAR_LSU LSU_reads
: > mapping_summary_SSU.txt
: > mapping_summary_LSU.txt

while IFS= read -r sample; do
  [[ -z "${sample}" ]] && continue
  filename="$(basename "${sample}")"
  echo ">>> Sample: ${filename}"

  # --- SSU mapping ---
  STAR \
    --genomeDir "${SSU_IDX}" \
    --readFilesIn "${sample}" \
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

  # SSU BAM path
  ssu_bam="SILVA_NR99_STAR_SSU/${filename}_Aligned.sortedByCoord.out.bam"

  # calculate the SSU mapping ratio
  ssu_mapped=$(samtools view -c -F 4   -F 0x100 -F 0x800 "${ssu_bam}")
  ssu_unmap=$(samtools view -c -f 4    -F 0x100 -F 0x800 "${ssu_bam}")
  printf "%s\t%s\t%s\n" "${filename}" "${ssu_mapped}" "${ssu_unmap}" >> mapping_summary_SSU.txt

  # extract the unmapped reads for LSU input
  lsu_fastq="LSU_reads/${filename}"
  samtools fastq -@ 4 -f 4 -F 0x100 -F 0x800 "${ssu_bam}" > "${lsu_fastq}"

  # --- LSU mapping ---
  STAR \
    --genomeDir "${LSU_IDX}" \
    --readFilesIn "${lsu_fastq}" \
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

  lsu_bam="SILVA_NR99_STAR_LSU/${filename}_Aligned.sortedByCoord.out.bam"
  lsu_mapped=$(samtools view -c -F 4   -F 0x100 -F 0x800 "${lsu_bam}")
  lsu_unmap=$(samtools view -c -f 4    -F 0x100 -F 0x800 "${lsu_bam}")
  printf "%s\t%s\t%s\n" "${filename}" "${lsu_mapped}" "${lsu_unmap}" >> mapping_summary_LSU.txt

done < "${LIST_20K}"

echo "Mapping finished at: $(date)"

# ========= 5) BLASTN verify  =========
cd "${SUBSET_DIR}"

python /mnt/data3/yiyonghao/MicroRNA/code_upload/4.rRNAcalculated/rRNA_blast.py \
  --input-dir ./SILVA_NR99_STAR_LSU \
  --db /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/blastn_core_nt/core_nt \
  --outdir /mnt/4T_SSD/yiyonghao/blastn_out \
  --samtools-threads 64 --blastn-threads 64 \
  --tmpdir /mnt/4T_SSD/yiyonghao/tmp


python /mnt/data3/yiyonghao/MicroRNA/code_upload/4.rRNAcalculated/rRNA_blastn_mapped.py
echo "All done: $(date)"
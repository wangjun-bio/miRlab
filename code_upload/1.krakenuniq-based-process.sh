# PS: this script is used to generated abundance files based on Krakenuniq
# modification : 2025-08-11
date

echo "Start making STAR index for human genome"
cd /mnt/data3/yiyonghao/genome-index/human-genome
STAR --runThreadN 64 \
     --runMode genomeGenerate \
     --genomeDir ./STAR_index \
     --genomeFastaFiles ./GCF_000001405.26_GRCh38_genomic.fna \
     --sjdbGTFfile ./GCF_000001405.40_GRCh38.p14_genomic.gtf \
     --sjdbOverhang 99


#### Trimmed 36 nt reads ####
# NC paper
echo "Start trimming 36 nt reads of NC paper"
cd /mnt/data3/yiyonghao/NC_paper_rawdata/upload
mkdir -p ./trimmed_36nt
mkdir -p ./decompressed
while IFS= read -r line
do 
    pigz -dc -p 64 ./${line}.fastq.gz > ./decompressed/${line}.fastq
done < NC_sample.txt

cat NC_sample.txt | parallel -j 40 trimmomatic SE ./decompressed/{}.fastq trimmed_36nt/{}.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36

# remove decompressed files
#rm -rf ./decompressed

# PN 
echo "Start trimming 36 nt reads of PN"
cd /mnt/data3/yiyonghao/NC_paper_rawdata/PN
find . -type f -name "*.gz" -exec basename {} .fastq.gz \; > sample.txt
mkdir -p ./trimmed_36nt
mkdir -p ./decompressed
while IFS= read -r line
do
    pigz -dc -p 64 ./${line}_cut_trim_R1.fastq.gz > ./decompressed/${line}.fastq
done < sample.txt
cat sample.txt | parallel -j 40 trimmomatic SE ./decompressed/{}.fastq trimmed_36nt/{}.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36
# remove decompressed files
rm -rf ./decompressed

# Brain
echo "Start trimming 36 nt reads of Brain"
cd /mnt/data3/yiyonghao/NC_paper_rawdata/Brain
find . -type f -name "*.gz" -exec basename {} .fastq.gz \; > sample.txt
mkdir -p ./trimmed_36nt
mkdir -p ./decompressed
while IFS= read -r line
do
    pigz -dc -p 64 ./${line}_cut_trim_R1.fastq.gz > ./decompressed/${line}.fastq
done < sample.txt
cat sample.txt | parallel -j 40 trimmomatic SE ./decompressed/{}.fastq trimmed_36nt/{}.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36
# remove decompressed files
rm -rf ./decompressed

# H20
echo "Start trimming 36 nt reads of H20"
cd /mnt/data3/yiyonghao/NC_paper_rawdata/h20_all
find . -type f -name "*.gz" -exec basename {} .fastq.gz \; > sample.txt
mkdir -p ./trimmed_36nt
mkdir -p ./decompressed
while IFS= read -r line
do
    pigz -dc -p 64 ./${line}_cut_trim_R1.fastq.gz > ./decompressed/${line}.fastq
done < sample.txt
cat H20_sample.txt | parallel -j 40 trimmomatic SE ./decompressed/{}.fastq trimmed_36nt/{}.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:36
# remove decompressed files
rm -rf ./decompressed



### STAR mapping to hg38 genome ####
#NC paper
cd /mnt/data3/yiyonghao/NC_paper_rawdata/upload/
mkdir -p ./star_mapping
echo ----------------------------------------------------
echo "Start mapping 36nt reads"
date
while IFS= read -r i; do
    STAR --runThreadN 60 \
        --genomeDir /mnt/data3/yiyonghao/genome-index/human-genome/STAR_index \
        --readFilesIn trimmed_36nt/${i}.fastq \
        --outSAMunmapped Within \
        --outFileNamePrefix star_mapping/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard
done < 'NC_sample.txt'
echo completed mapping 36nt reads of NC paper
date

# PN
cd /mnt/data3/yiyonghao/NC_paper_rawdata/PN/
mkdir -p ./star_mapping
echo ----------------------------------------------------
echo "Start mapping 36nt reads"
date
while IFS= read -r i; do
    STAR --runThreadN 60 \
        --genomeDir /mnt/data3/yiyonghao/genome-index/human-genome/STAR_index \
        --readFilesIn trimmed_36nt/${i}.fastq \
        --outSAMunmapped Within \
        --outFileNamePrefix star_mapping/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard
done < 'sample.txt'
echo completed mapping 36nt reads of PN
date

# Brain
cd /mnt/data3/yiyonghao/NC_paper_rawdata/Brain/
mkdir -p ./star_mapping
echo ----------------------------------------------------
echo "Start mapping 36nt reads"
date
while IFS= read -r i; do
    STAR --runThreadN 60 \
        --genomeDir /mnt/data3/yiyonghao/genome-index/human-genome/STAR_index \
        --readFilesIn trimmed_36nt/${i}.fastq \
        --outSAMunmapped Within \
        --outFileNamePrefix star_mapping/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard
done < 'sample.txt'
echo completed mapping 36nt reads of Brain
date

# H20
cd /mnt/data3/yiyonghao/NC_paper_rawdata/h20_all/
mkdir -p ./star_mapping
echo ----------------------------------------------------
echo "Start mapping 36nt reads"
date
while IFS= read -r i; do
    STAR --runThreadN 60 \
        --genomeDir /mnt/data3/yiyonghao/genome-index/human-genome/STAR_index \
        --readFilesIn trimmed_36nt/${i}.fastq \
        --outSAMunmapped Within \
        --outFileNamePrefix star_mapping/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard
done < 'sample.txt'
echo completed mapping 36nt reads of H20
date



## Extract unmapped reads from STAR output ##
# NC paper
cd /mnt/data3/yiyonghao/NC_paper_rawdata/upload/star_mapping/
mkdir -p ./unmapped_reads
echo ----------------------------------------------------
echo "Start extracting unmapped reads of NC paper"
date
while IFS= read -r i; do
    samtools view -@ 64 -b -f 4 ${i}_Aligned.sortedByCoord.out.bam > unmapped_reads/${i}_unmapped.bam
    samtools fastq -@ 64 unmapped_reads/${i}_unmapped.bam > unmapped_reads/${i}_unmapped.fastq
done < '../sample.txt'
rm -rf ./unmapped_reads/*.bam
echo completed extracting unmapped reads of NC paper
date

# PN
cd /mnt/data3/yiyonghao/NC_paper_rawdata/PN/star_mapping
mkdir -p ./unmapped_reads
echo ----------------------------------------------------
echo "Start extracting unmapped reads of PN"
date
while IFS= read -r i; do
    samtools view -@ 64 -b -f 4 ${i}_Aligned.sortedByCoord.out.bam > unmapped_reads/${i}_unmapped.bam
    samtools fastq -@ 64 unmapped_reads/${i}_unmapped.bam > unmapped_reads/${i}_unmapped.fastq
done < '../sample.txt'
rm -rf ./unmapped_reads/*.bam
echo completed extracting unmapped reads of PN
date

# Brain
cd /mnt/data3/yiyonghao/NC_paper_rawdata/Brain/star_mapping
mkdir -p ./unmapped_reads
echo ----------------------------------------------------
echo "Start extracting unmapped reads of Brain"
date
while IFS= read -r i; do
    samtools view -@ 64 -b -f 4 ${i}_Aligned.sortedByCoord.out.bam > unmapped_reads/${i}_unmapped.bam
    samtools fastq -@ 64 unmapped_reads/${i}_unmapped.bam > unmapped_reads/${i}_unmapped.fastq
done < '../sample.txt'
rm -rf ./unmapped_reads/*.bam
echo completed extracting unmapped reads of Brain
date


# H20
cd /mnt/data3/yiyonghao/NC_paper_rawdata/h20_all/star_mapping
mkdir -p ./unmapped_reads
echo ----------------------------------------------------
echo "Start extracting unmapped reads of H20"
date
while IFS= read -r i; do
    samtools view -@ 64 -b -f 4 ${i}_Aligned.sortedByCoord.out.bam > unmapped_reads/${i}_unmapped.bam
    samtools fastq -@ 64 unmapped_reads/${i}_unmapped.bam > unmapped_reads/${i}_unmapped.fastq
done < '../sample.txt'
rm -rf ./unmapped_reads/*.bam
echo completed extracting unmapped reads of H20
date




## Krakenuniq mapping ##
# NC paper
cd /mnt/data3/yiyonghao/NC_paper_rawdata/upload/
mkdir -p ./krakenuniq_output

/mnt/data3/wangjun/software/krakenuniq-1.0.4/krakenuniq --db /mnt/data3/wangjun/cmRNA_20241107/ref/krakenuniq_db_20240122/krakenuniq_standard --preload --threads 64

echo ----------------------------------------------------
echo "Start Krakenuniq classification for NC paper"
cat sample.txt | parallel -j 64 '/mnt/data3/wangjun/software/krakenuniq-1.0.4/krakenuniq --db /mnt/data3/wangjun/cmRNA_20241107/ref/krakenuniq_db_20240122/krakenuniq_standard \
    --report-file ./krakenuniq_output/{}_reportfile.tsv ./star_mapping/unmapped_reads/{}_unmapped.fastq > ./krakenuniq_output/{}_readclassification.tsv'

# PN
cd /mnt/data3/yiyonghao/NC_paper_rawdata/PN/
mkdir -p ./krakenuniq_output
echo ----------------------------------------------------
echo "Start Krakenuniq classification for PN"
cat sample.txt | parallel -j 64 '/mnt/data3/wangjun/software/krakenuniq-1.0.4/krakenuniq --db /mnt/data3/wangjun/cmRNA_20241107/ref/krakenuniq_db_20240122/krakenuniq_standard \
    --report-file ./krakenuniq_output/{}_reportfile.tsv ./star_mapping/unmapped_reads/{}_unmapped.fastq > ./krakenuniq_output/{}_readclassification.tsv'

# Brain
cd /mnt/data3/yiyonghao/NC_paper_rawdata/Brain/
mkdir -p ./krakenuniq_output
echo ----------------------------------------------------
echo "Start Krakenuniq classification for Brain"
cat sample.txt | parallel -j 64 '/mnt/data3/wangjun/software/krakenuniq-1.0.4/krakenuniq --db /mnt/data3/wangjun/cmRNA_20241107/ref/krakenuniq_db_20240122/krakenuniq_standard \
    --report-file ./krakenuniq_output/{}_reportfile.tsv ./star_mapping/unmapped_reads/{}_unmapped.fastq > ./krakenuniq_output/{}_readclassification.tsv'


# PN
cd /mnt/data3/yiyonghao/NC_paper_rawdata/h20_all/
mkdir -p ./krakenuniq_output
echo ----------------------------------------------------
echo "Start Krakenuniq classification for H20"
cat sample.txt | parallel -j 64 '/mnt/data3/wangjun/software/krakenuniq-1.0.4/krakenuniq --db /mnt/data3/wangjun/cmRNA_20241107/ref/krakenuniq_db_20240122/krakenuniq_standard \
    --report-file ./krakenuniq_output/{}_reportfile.tsv ./star_mapping/unmapped_reads/{}_unmapped.fastq > ./krakenuniq_output/{}_readclassification.tsv'

echo "All processes completed successfully."
date

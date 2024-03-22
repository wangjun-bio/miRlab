# cfMeDIP-seq pipline
# move raw fastq file to current folder
# 20240322

echo ------
echo note: prepare dir and file...
date
mkdir rawdata result_fastqc result_cutadapter result_trimmomatic mappeddata_bwa
cp barcode.3.fasta ./
cp split_barcode.py ./rawdata/
cp dedup.py ./mappeddata_bwa/
cp unique_mapped_filter.py ./mappeddata_bwa/


echo ------
echo note: creat the index of sample list...
date
ls *_R1.fastq | cut -d '_' -f 1 > sample_list.txt


echo ------
echo note: move the fastq into the rawdata
date
mv *.fastq rawdata

echo ------
echo note: begin fastqc
date
fastqc -q -t 40 -o ./result_fastqc/ ./rawdata/*.fastq


echo ------
echo note: cut 5 bar begin
date
cd rawdata
for i in `cat ../sample_list.txt`
do
	python split_barcode.py ${i}_R1.fastq ${i}_R2.fastq
done
cd ..
echo note: cut 5 bar done
date


echo ------
echo note: begin cut adapter...
date
for i in `cat sample_list.txt`
do
	cutadapt -j 40 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ./result_cutadapter/${i}_cut5bar_cutadapt_R1.fastq -p ./result_cutadapter/${i}_cut5bar_cutadapt_R2.fastq ./rawdata/${i}_split_barcode_R1.fastq ./rawdata/${i}_split_barcode_R2.fastq
done
echo note: cut adapter done
date


echo ------
echo note: begin cut 3 bar...
date
for i in `cat sample_list.txt`
do
	cutadapt -j 40 -e 0 --no-indels -a file:barcode.3.fasta -A file:barcode.3.fasta -o ./result_cutadapter/${i}_cut5bar_cutadapt_cut3bar_R1.fastq -p ./result_cutadapter/${i}_cut5bar_cutadapt_cut3bar_R2.fastq ./result_cutadapter/${i}_cut5bar_cutadapt_R1.fastq ./result_cutadapter/${i}_cut5bar_cutadapt_R2.fastq
done
echo note: cut 3 bar done
date


echo ------
echo note: begin trim...
date
for i in `cat sample_list.txt`
do
	trimmomatic PE -threads 40 ./result_cutadapter/${i}_cut5bar_cutadapt_cut3bar_R1.fastq ./result_cutadapter/${i}_cut5bar_cutadapt_cut3bar_R2.fastq -baseout ./result_trimmomatic/${i}_cut5bar_cutadapt_cut3bar_trim.fastq LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20
done
echo note: trim done
date


echo ------
echo note: begin mapping...
for i in `cat sample_list.txt`
do
	 bwa mem -t 40 -M /home/Ref/Hg19/hg19.fa ./result_trimmomatic/${i}_cut5bar_cutadapt_cut3bar_trim_1P.fastq ./result_trimmomatic/${i}_cut5bar_cutadapt_cut3bar_trim_2P.fastq > ./mappeddata_bwa/${i}.sam
done
echo note: mapping done
date


echo ------
echo note: get unique mapped reads...
date
cd mappeddata_bwa
for i in `cat ../sample_list.txt`
do
	nohup python unique_mapped_filter.py ${i}.sam
done
cd ..
echo note: uniqmad done
date


echo ------
echo note: remove the duplication...
date
cd mappeddata_bwa
for i in `cat ../sample_list.txt`
do
        nohup python dedup.py ${i}_uniqmap.sam
done
cd ..
echo note: dedup done
date


echo ------
echo note: sam to bam, sort and build index...
date
for i in `cat sample_list.txt`
do
        samtools view -@ 40 -S ./mappeddata_bwa/${i}_uniqmap_dedup.sam -b > ./mappeddata_bwa/${i}_uniqmap_dedup.bam  
        samtools sort -@ 40 ./mappeddata_bwa/${i}_uniqmap_dedup.bam -o ./mappeddata_bwa/${i}_uniqmap_dedup_sorted.bam   
        samtools index -@ 40 ./mappeddata_bwa/${i}_uniqmap_dedup_sorted.bam   
        
        samtools view -@ 40 -Sb ./mappeddata_bwa/${i}.sam > ./mappeddata_bwa/${i}.bam
        samtools flagstat -@ 40 ./mappeddata_bwa/${i}.bam > ./mappeddata_bwa/${i}_flagstat.txt
done
echo note: bam done
date



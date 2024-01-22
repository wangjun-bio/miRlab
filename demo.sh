#!/bin/bash 
echo start================================================================================================
date

path_raw_read=./demo
path_raw_cutadapt=./raw_cutadapt
path_raw_trimmomatic=./raw_trimmomatic

#cutadapt
echo start cutadapt=========================================================================================
date
for i in $(cat $path_raw_read/samplelist.txt); 
do 
    echo  cut +++++++++++++ ${i}-start
    cutadapt -j 20 -a AGATCGGAAAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -O 5 -e 0 -o $path_raw_cutadapt/${i}_cutadapt_1.fastq.gz -p $path_raw_cutadapt/${i}_cutadapt_2.fastq.gz $path_raw_read/${i}_1.fq.gz $path_raw_read/${i}_2.fq.gz 
    echo   ${i}_cutadapt_1.fastq ========finsh cutadapt
done
#trimmomatic
echo start trimmomatic=====================================================================================
date
for i in $(cat $path_raw_read/samplelist.txt); 
do 
    echo trim+++++++++++ ${i}-start
    trimmomatic PE -threads 20 $path_raw_cutadapt/${i}_cutadapt_1.fastq.gz $path_raw_cutadapt/${i}_cutadapt_2.fastq.gz -baseout $path_raw_trimmomatic/${i}_cutam.fastq.gz LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 AVGQUAL:20 MINLEN:20;
    echo ${i}_cutam.fastq =========finsh trimmomatic
done


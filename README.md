# Introduction
Scripts in this folder are collection of bioinformatics analysis piplines to process cfMeDIP-seq data. These piplines take paired-end FASTQ files as input and generate output files which can be used for downstream analysis and the cfDNA fragment size calculation. 

# upstream_pipline
Scripts in this part could be used to generate bam file from raw fastq file.<br>raw fastq files and script "cfMeDIP-seq_upstream.sh"， "barcode.3.fasta", "dedup.py","split_barcode.py","unique_mapped_filter.py" shoul be placed in the same folder. <br> Through runing the cfMeDIP-seq_upstream script will result the bam file for downstream analysis.<br>test_R1.fastq and test_R2.fastq could be used as demo data. 

# downstream_pipline
Scripts in this folder could be used to calculate cfDNA fragment size in different group of individuals.<br>01.frags_generation.R, this script use the bame file generated from unstream_pipline as input, and then generate "frags.Rdata" file.<br>02.IP_Input_fragment_plot.R, this script use the "test_group_info.csv" and bin size, such as "10kb" as input.<br>
03.DMR_generation.R, this script use the "test_group.csv" and bin size, such as "10kb" as input to identify DMR, and calculate fragment size distribution with DMR region.<br>04.DMR_plot.R, this script use the Rdata file generated by 02.IP_Input_fragment_plot.R script as input, at the same time, padj_val such as "0.05" and log2fold_cal such as "0.5" are also used as input.
# miRlab
## Data Availability

### Raw sequencing data
The raw FASTQ files have been deposited in the Gene Expression Omnibus (GEO) under accession number **GSEXXXXX** (will be released upon publication).

### Processed data (expression matrices, mutation calls, microbial abundances)
Due to GitHub storage limitations, all processed data files (including RPKM matrices, count tables, mutation presence matrices, DE results, and featureCounts summaries) are available upon request and shoule be place in the following directories relative to the project root:

- `data/df_rpkm_rpm.rds`
- `data/df_counts.rds`
- `data/df_rpkm_300W_allRNA_microbe_variant.csv`
- `data/DE_results/` (all DE CSV files)
- `data/annovar/csv_file/` (ANNOVAR output CSV files)
- `data/krakenuniq/20250707/` (KrakenUniq report files)
- `results/featureCounts/summaries/` (mRNA/lncRNA/snRNA/snoRNA assigned counts)

### Reference genome and annotation files
The following reference files can be obtained from public sources 
- **GENCODE v38 GTF** (Coding_gene_annotation.gtf, lncRNA annotation, etc.):  
- **Bowtie2 index for hg38**:  
  `bowtie2-build hg38.fa hg38`
- **KrakenUniq standard database**:  

### Sample lists and RNA name files
These small text files are already included in the repository under `data/`.

### Running the analysis from scratch
All analysis steps can be fully reproduced from the raw FASTQ files using the scripts in `scripts/`. Please see `README.md` for the complete workflow.
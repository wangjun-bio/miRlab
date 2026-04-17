# miRlab
# PAH cfRNA multi-omics analysis

This repository contains the complete analysis pipeline for the paper "".

## Data availability
Raw sequencing data are available at GEO under accession GSEXXXXX.  
Processed expression matrices and sample annotations can be found in `data/` (not included due to size; contact author for access).

## Dependencies
- R 4.2+ with packages: DESeq2, ggplot2, pheatmap, survival, maftools, etc.
- Python 3.9+ with: pandas, numpy, scikit-learn, pydeseq2, lifelines, etc.
- Conda environment: `conda env create -f environment.yml`

## Usage
1. Clone the repo
2. Prepare input data (see `data/README.md`)
3. Run analysis in order:
   - `bash/01_trim_align.sh`
   - `R/01_QC_visualization.R`
   - ...

## License
MIT
#!/usr/bin/env bash
IFS=$'\n\t'

# check tools 
for tool in awk shuf seqtk; do
  command -v $tool >/dev/null 2>&1 || { echo >&2 "Error: '$tool' not found in PATH."; exit 1; }
done

# Paths
dirs=(
  "/mnt/data3/yiyonghao/NC_paper_rawdata/upload"
  "/mnt/data3/yiyonghao/NC_paper_rawdata/Brain"
  "/mnt/data3/yiyonghao/NC_paper_rawdata/PN/"
)

# index files
index_file="sample.txt"

# Traverse the all path
for d in "${dirs[@]}"; do
  echo "=== Processing directory: $d ==="
  (
    cd "$d"
    echo "Current working dir: $(pwd)"
    mkdir -p ./blastn/info
    mkdir -p ./blastn/index
    mkdir -p ./blastn/fastq
    mkdir -p ./blastn/results
    # check the index existed 
    if [[ ! -f "$index_file" ]]; then
      echo "Warning: Index file '$index_file' not found in $d; skip this directory." >&2
      continue
    fi

    # randomly select reads
    while IFS= read -r sample; do
      [[ -z "$sample" ]] && continue
      tsv="./krakenuniq_output/${sample}_readclassification.tsv"
      fq1="./decompressed/${sample}.fastq"
      if [[ -f "$fq1" ]]; then
        fastq="$fq1"
      else
        echo "Warning: neither '$fq1' nor found; skipping $sample" >&2
        continue
      fi

      out_info="./blastn/info/${sample}_random100_C.tsv"
      id_list="./blastn/index/${sample}_random100_ids.txt"
      out_fastq="./blastn/fastq/${sample}_random100.fastq"

      echo "--- Sample: $sample ---"
      # 1) filter C row and randomly select 100 reads
      awk '$1=="C"' "$tsv" | shuf -n 100 > "$out_info"
      echo "  -> saved classification lines to $out_info"

      # 2) extract the id infomation
      awk '{print $2}' "$out_info" > "$id_list"
      echo "  -> extracted IDs to $id_list"

      # 3) extract reads from FASTQ
      seqtk subseq "$fastq" "$id_list" > "$out_fastq"
      echo "  -> saved sequences to $out_fastq"
      echo "  -> Processed $sample successfully."
      rm "$id_list"
    done < "$index_file"
  )
 ## Merge the all FASTQ
  echo "Merging FASTQ files in $d/blastn/fastq"
  merged_fastq="$d/blastn/fastq/merged_random100.fastq"
  merged_fasta="$d/blastn/fastq/merged_random100.fasta"
  if ls "$d/blastn/fastq/"*.fastq >/dev/null 2>&1; then
    cat "$d/blastn/fastq/"*.fastq > "$merged_fastq"
    seqtk seq -A $merged_fastq > $merged_fasta
    echo "  -> Merged FASTQ saved to $merged_fastq"
  else
    echo "Warning: No FASTQ files found in $d/blastn/fastq; skipping merge." >&2
  fi  
  
  # execute Blastn
  echo "Running blastn on merged FASTA: $merged_fasta"
  blastn_output="$d/blastn/results/merged_random100_blastn.out"
  blastn -query "$merged_fasta" -db /mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/blastn_core_nt/core_nt \
    -out "$blastn_output" -outfmt 6 -num_threads 60
  echo "  -> Blastn results saved to $blastn_output"
  echo "=== Done directory: $d ===\n"
done

echo "All directories processed."

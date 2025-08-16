rm -rf ./length
mkdir -p length
rm -rf ./sample.txt
find . -name "*.fastq" -exec basename {} \; | sed 's/\.fastq$//' > sample.txt
while IFS= read -r line;
do
	echo "${line} is processing"
	python calculateLength.py -i ./${line}.fastq -o ./length/${line}.tsv -t 60
done < "sample.txt"

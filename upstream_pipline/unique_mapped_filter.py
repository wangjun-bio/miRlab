
# keep the uniquely mapped reads in the bwa mem produced sam file
import os
import linecache
import re
import sys
count=0
# input file and define the output file name and path for save
sam=sys.argv[1]
name=sam.split('.')[0].strip()
file_path=os.getcwd()
output=open(file_path+'/'+name+'_uniqmap.sam','w')
sam_1=open(sam,'r')
# read every line of sam file starting from 1st line and omit the information stored in the starting region of the sam file, and calculate the total reads number
for line in sam_1:
	if not line.startswith('@'):
		count+=1
#print (str(count)+'reads will be processed')
sam_2=open(sam)
# define the flag dictionary, and only save the flag=99,147,83,163, which are reads that have uniq map
d_multiple={}
d_uniq={}
uniq_map_count=0
multiple_map_count=0
head={}
# read every line of sam file, and save read_ID, read_MB,read_flag,read_mapQ,read_chr,read_start,read_end,insert_size
for line in sam_2:
	if line.startswith('@'):
		head_info=line
		head[head_info]=head_info
	elif not line.startswith('@'):
		mapping_info=line
		read_ID=mapping_info.strip().split('\t')[0].split(':B_')[0]
		multiple_map_marker=mapping_info.strip().split('\t')[-1]
		if re.match('.A:Z:',multiple_map_marker):
			d_multiple[read_ID]=multiple_map_marker
sam_3=open(sam)
for line in sam_3:
	if not line.startswith('@'):
		mapping_info=line
		read_ID=mapping_info.strip().split('\t')[0].split(':B_')[0]
		read_flag=str(mapping_info.strip().split('\t')[1])
		read_chr=mapping_info.strip().split('\t')[2]
		read_start=mapping_info.strip().split('\t')[3]
		read_end=mapping_info.strip().split('\t')[7]
		insert_size=mapping_info.strip().split('\t')[8]
		if read_ID not in d_multiple.keys():
			reads_key=tuple([str(read_ID),str(read_flag),str(read_chr),str(read_start),str(read_end),str(insert_size)])
			d_uniq[reads_key]=mapping_info
			uniq_map_count+=1
		elif read_ID in d_multiple.keys():
			multiple_map_count+=1
print('writting',name,'to output files...')
# save thr unduplicated reads into a new file
for keys in head.keys():
	output.write(head[keys])
for keys in d_uniq.keys():
	output.write(d_uniq[keys])
print(str(uniq_map_count) + ' uniq mapped reads have been found in',name)
print(str(multiple_map_count)+' multiple mapped reads have been found in',name)
output.close()

# use molecular barcode to remove duplicaiton in sam file
# calculate the total reads of input sam file
import os
import linecache
import re
import sys
count=0
# input file and define the output file name and path for save
sam=sys.argv[1]
name=sam.split('.')[0].strip()
file_path=os.getcwd()
output=open(file_path+'/'+name+'_dedup.sam','w')
sam_1=open(sam,'r')
# read every line of sam file starting from 1st line and omit the information stored in the starting region of the sam file, and calculate the total reads number
for line in sam_1:
	if not line.startswith('@'):
		count+=1
#print (str(count)+'reads will be processed')
sam_2=open(sam)
# define the flag dictionary, and only save the flag=99,147,83,163, which are reads that have properly paired map
d_99={}
d_147={}
d_83={}
d_163={}
head={}
properly_paired_flag=['99','147','83','163']
properly_paired_reads=0
dup_count=0
# read every line of sam file, and save read_ID, read_MB,read_flag,read_mapQ,read_chr,read_start,read_end,insert_size
for line in sam_2:
	if line.startswith('@'):
		head_info=line
		head[head_info]=head_info
	elif not line.startswith('@'):
		mapping_info=line
		read_ID=mapping_info.strip().split('\t')[0].split(':B_')[0]
		read_MB=mapping_info.strip().split('\t')[0].split(':')[7]
		read_flag=str(mapping_info.strip().split('\t')[1])
		read_mapQ=int(mapping_info.strip().split('\t')[4])
		read_chr=mapping_info.strip().split('\t')[2]
		read_start=mapping_info.strip().split('\t')[3]
		read_end=mapping_info.strip().split('\t')[7]
		insert_size=mapping_info.strip().split('\t')[8]
		# first process the flag=99 reads. for these reads, using molecular-barcode,chr,start,insert-size for looking for duplication reads  
		if read_flag == '99':
			dup_info=tuple([str(read_MB),str(read_chr),str(read_start),str(insert_size)])
			properly_paired_reads+=1
			if dup_info in d_99.keys():
				dup_count+=1
				#if read_mapQ == int(d_99[dup_info].strip().split('\t')[4]):
				#	print ('Remove:'+read_ID+'---- same mapQ found ----dupfrom:'+d_99[dup_info].strip().split('\t')[0].split(':B_')[0])
				#elif read_mapQ < int(d_99[dup_info].strip().split('\t')[4]):
				#	print('Remove:' +read_ID+'---- lower mapQ found ----dupfrom:'+d_99[dup_info].strip().split('\t')[0].split(':B_')[0])
				## save the higher mapQ reads if duplication found
				if read_mapQ > int(d_99[dup_info].strip().split('\t')[4]):
				#	print('Remove:'+d_99[dup_info].strip().split('\t')[0].split(':B_')[0]+'---- higher mapQ found----dupfrom:'+read_ID)
					d_99[dup_info]=mapping_info
			elif dup_info not in d_99.keys():
				d_99[dup_info]=mapping_info
		# process the flag=147 reads, for these reads, using molecular-barcode,chr,end,insert-size for looking for duplication reads
		elif read_flag == '147':
			dup_info=tuple([str(read_MB),str(read_chr),str(read_end),str(insert_size)])
			properly_paired_reads+=1
			if dup_info in d_147.keys():
				dup_count+=1
				#if read_mapQ == int(d_147[dup_info].strip().split('\t')[4]):
				#	print ('Remove:'+read_ID+'---- same mapQ found ----dupfrom:'+d_147[dup_info].strip().split('\t')[0].split(':B_')[0])
				#elif read_mapQ < int(d_147[dup_info].strip().split('\t')[4]):
				#	print('Remove:' +read_ID+'---- lower mapQ found ----dupfrom:'+d_147[dup_info].strip().split('\t')[0].split(':B_')[0])
				if read_mapQ > int(d_147[dup_info].strip().split('\t')[4]):
				#	print('Remove:'+d_147[dup_info].strip().split('\t')[0].split(':B_')[0]+'---- higher mapQ found----dupfrom:'+read_ID)
					d_147[dup_info]=mapping_info
			elif dup_info not in d_147.keys():
				d_147[dup_info]=mapping_info
		# process the flag=163 reads,similar to flag=99
		elif read_flag == '163':
			dup_info=tuple([str(read_MB),str(read_chr),str(read_start),str(insert_size)])
			properly_paired_reads+=1
			if dup_info in d_163.keys():
				dup_count+=1
				#if read_mapQ == int(d_163[dup_info].strip().split('\t')[4]):
				#	print ('Remove:'+read_ID+'---- same mapQ found ----dupfrom:'+d_163[dup_info].strip().split('\t')[0].split(':B_')[0])
				#elif read_mapQ < int(d_163[dup_info].strip().split('\t')[4]):
				#	print('Remove:' +read_ID+'---- lower mapQ found ----dupfrom:'+d_163[dup_info].strip().split('\t')[0].split(':B_')[0])
				if read_mapQ > int(d_163[dup_info].strip().split('\t')[4]):
				#	print('Remove:'+d_163[dup_info].strip().split('\t')[0].split(':B_')[0]+'---- higher mapQ found----dupfrom:'+read_ID)
					d_163[dup_info]=mapping_info
			elif dup_info not in d_163.keys():
				d_163[dup_info]=mapping_info
		# process the flag=83 reads, similar to flag=147		
		elif read_flag == '83':
			dup_info=tuple([str(read_MB),str(read_chr),str(read_end),str(insert_size)])
			properly_paired_reads+=1
			if dup_info in d_83.keys():
				dup_count+=1 
				#if read_mapQ == int(d_83[dup_info].strip().split('\t')[4]):
				#	print ('Remove:'+read_ID+'---- same mapQ found ----dupfrom:'+d_83[dup_info].strip().split('\t')[0].split(':B_')[0])
				#elif read_mapQ < int(d_83[dup_info].strip().split('\t')[4]):
				#	print('Remove:' +read_ID+'---- lower mapQ found ----dupfrom:'+d_83[dup_info].strip().split('\t')[0].split(':B_')[0])
				if read_mapQ > int(d_83[dup_info].strip().split('\t')[4]):
				#	print('Remove:'+d_83[dup_info].strip().split('\t')[0].split(':B_')[0]+'---- higher mapQ found----dupfrom:'+read_ID)
					d_83[dup_info]=mapping_info
			elif dup_info not in d_83.keys():
				d_83[dup_info]=mapping_info
print('writting to output files...')
# save thr unduplicated reads into a new file
for keys in head.keys():
	output.write(head[keys])
for keys in d_99.keys():
	output.write(d_99[keys])
for keys in d_147.keys():
	output.write(d_147[keys])
for keys in d_163.keys():
	output.write(d_163[keys])
for keys in d_83.keys():
	output.write(d_83[keys])
print(str(count) + ' PE reads have been processed')
print(str(properly_paired_reads)+' PE reads are properly paired among total reads processed')
print(str(dup_count)+' PE reads are duplication among  properly paired reads')
print(str(properly_paired_reads-dup_count)+" PE reads are unduplicated properly paired reads")
output.close()

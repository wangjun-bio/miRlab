import re
import sys
import linecache
import os
# 输入R1与R2文件
R1=sys.argv[1]
R2=sys.argv[2]
#用字典定义八个不同的分子标签，并允许分子标签中有一个错配
MB={}
MB['B_1']=['CATGCCTA','.ATGCCTA','C.TGCCTA','CA.GCCTA','CAT.CCTA','CATG.CTA','CATGC.TA','CATGCC.A','CATGCCT.']
MB['B_2']=['AGCGTAGC','.GCGTAGC','A.CGTAGC','AG.GTAGC','AGC.TAGC','AGCG.AGC','AGCGT.GC','AGCGTA.C','AGCGTAG.']
MB['B_3']=['GCTCAGGA','.CTCAGGA','G.TCAGGA','GC.CAGGA','GCT.AGGA','GCTC.GGA','GCTCA.GA','GCTCAG.A','GCTCAGG.']
MB['B_4']=['TTCTGCCT','.TCTGCCT','T.CTGCCT','TT.TGCCT','TTC.GCCT','TTCT.CCT','TTCTG.CT','TTCTGC.T','TTCTGCC.']
MB['B_5']=['CGCTCTTC','.GCTCTTC','C.CTCTTC','CG.TCTTC','CGC.CTTC','CGCT.TTC','CGCTC.TC','CGCTCT.C','CGCTCTT.']
MB['B_6']=['TATGCGTT','.ATGCGTT','T.TGCGTT','TA.GCGTT','TAT.CGTT','TATG.GTT','TATGC.TT','TATGCG.T','TATGCGT.']
MB['B_7']=['GGACCGAT','.GACCGAT','G.ACCGAT','GG.CCGAT','GGA.CGAT','GGAC.GAT','GGACC.AT','GGACCG.T','GGACCGA.']
MB['B_8']=['CCAGGTTC','.CAGGTTC','C.AGGTTC','CC.GGTTC','CCA.GTTC','CCAG.TTC','CCAGG.TC','CCAGGT.C','CCAGGTT.']
#定义输出文件的名称，并将输出文件存放在当前文件夹
R1_name=str(R1.split('_')[0]+'_split_barcode_R1.fastq')
R2_name=str(R2.split('_')[0]+'_split_barcode_R2.fastq')
file_path=os.getcwd()
R1_output=open(file_path+'/'+R1_name,'a')
R2_output=open(file_path+'/'+R2_name,'a')
#首先打开R1
R1=open(R1,'r')
i=1
for line in R1:
	# 从第1行开始，依次读取R1.fq的每一行内容，并标记他们的类型
	if i%4==1:
		R1_seq_ID=line.split(' ')[0].strip('\n')
	elif i%4==2:
		R1_seq=line.strip('\n')
		# 提取R1最左侧的8bp作为分子标签，剩余的9-150个bp作为真正的测序序列
		R1_MB=R1_seq.strip('\n')[0:8]
		R1_true_seq=R1_seq.strip('\n')[9:151]
	elif i%4==3:
		R1_seq_plus=line.strip('\n')
	elif i%4==0:
		R1_MB_qua=line.strip('\n')[0:8]
		R1_true_seq_qua=line.strip('\n')[9:151]
		# 读取R2文件，此时i是4的整数倍，即此时刚读完一条reads的完整信息
		R2_seq_ID=linecache.getline(R2,i-3).split(' ')[0].strip('\n')
		R2_seq=linecache.getline(R2,i-2).strip('\n')
		R2_MB=R2_seq.strip('\n')[0:8]
		R2_true_seq=R2_seq.strip('\n')[9:151]
		R2_seq_plus=linecache.getline(R2,i-1).strip('\n')
		R2_MB_qua=linecache.getline(R2,i).strip('\n')[0:8]
		R2_true_seq_qua=linecache.getline(R2,i).strip('\n')[9:151]
		# 判断R1的分子标签是否存在于上述字典中，如果存在，记录是匹配到哪一个分子标签
		for keys in MB.keys():
			for item in MB[keys]:
				if re.match(item,R1_MB):
					#print (keys)
					R1_MB=str(keys)
					# 若R1成功匹配到某一个分子标签后，判断R2的分子标签是否也存在于上述字典中，记录是匹配到了哪一个分子标签
					for keys in MB.keys():
						for item in MB[keys]:
							if re.match(item,R2_MB):
								#print (keys)
								R2_MB=str(keys)
								#print (R1_MB+'-'+R2_MB)
			#若R1的分子标签和R2的分子标签同时都在上述字典中找到，则将该R1、R2的reads写入到输出文件中
		if R1_MB in MB.keys() and R2_MB in MB.keys():
			R1_seq_ID=str(R1_seq_ID+':'+str(R1_MB)+'+'+str(R2_MB))
			R2_seq_ID=str(R2_seq_ID+':'+str(R1_MB)+'+'+str(R2_MB))
			R1_output.write(R1_seq_ID+'\n'+R1_true_seq+'\n'+R1_seq_plus+'\n'+R1_true_seq_qua+'\n')
			R2_output.write(R2_seq_ID+'\n'+R2_true_seq+'\n'+R2_seq_plus+'\n'+R2_true_seq_qua+'\n')
	
		
	#运行过程中，每运行完1M条后，即打印出来		
	i+=1
	if i%4000000==0:
		print (str(int(i/4000000))+'M...reads processed')	
R1_output.close()
R2_output.close()

import re
import sys
import linecache
import os

R1=sys.argv[1]
R2=sys.argv[2]

MB={}
MB['B_1']=['CATGCCTA','.ATGCCTA','C.TGCCTA','CA.GCCTA','CAT.CCTA','CATG.CTA','CATGC.TA','CATGCC.A','CATGCCT.']
MB['B_2']=['AGCGTAGC','.GCGTAGC','A.CGTAGC','AG.GTAGC','AGC.TAGC','AGCG.AGC','AGCGT.GC','AGCGTA.C','AGCGTAG.']
MB['B_3']=['GCTCAGGA','.CTCAGGA','G.TCAGGA','GC.CAGGA','GCT.AGGA','GCTC.GGA','GCTCA.GA','GCTCAG.A','GCTCAGG.']
MB['B_4']=['TTCTGCCT','.TCTGCCT','T.CTGCCT','TT.TGCCT','TTC.GCCT','TTCT.CCT','TTCTG.CT','TTCTGC.T','TTCTGCC.']
MB['B_5']=['CGCTCTTC','.GCTCTTC','C.CTCTTC','CG.TCTTC','CGC.CTTC','CGCT.TTC','CGCTC.TC','CGCTCT.C','CGCTCTT.']
MB['B_6']=['TATGCGTT','.ATGCGTT','T.TGCGTT','TA.GCGTT','TAT.CGTT','TATG.GTT','TATGC.TT','TATGCG.T','TATGCGT.']
MB['B_7']=['GGACCGAT','.GACCGAT','G.ACCGAT','GG.CCGAT','GGA.CGAT','GGAC.GAT','GGACC.AT','GGACCG.T','GGACCGA.']
MB['B_8']=['CCAGGTTC','.CAGGTTC','C.AGGTTC','CC.GGTTC','CCA.GTTC','CCAG.TTC','CCAGG.TC','CCAGGT.C','CCAGGTT.']

R1_name=str(R1.split('_')[0]+'_split_barcode_R1.fastq')
R2_name=str(R2.split('_')[0]+'_split_barcode_R2.fastq')
file_path=os.getcwd()
R1_output=open(file_path+'/'+R1_name,'a')
R2_output=open(file_path+'/'+R2_name,'a')

R1=open(R1,'r')
i=1
for line in R1:
	
	if i%4==1:
		R1_seq_ID=line.split(' ')[0].strip('\n')
	elif i%4==2:
		R1_seq=line.strip('\n')
		
		R1_MB=R1_seq.strip('\n')[0:8]
		R1_true_seq=R1_seq.strip('\n')[9:151]
	elif i%4==3:
		R1_seq_plus=line.strip('\n')
	elif i%4==0:
		R1_MB_qua=line.strip('\n')[0:8]
		R1_true_seq_qua=line.strip('\n')[9:151]
		
		R2_seq_ID=linecache.getline(R2,i-3).split(' ')[0].strip('\n')
		R2_seq=linecache.getline(R2,i-2).strip('\n')
		R2_MB=R2_seq.strip('\n')[0:8]
		R2_true_seq=R2_seq.strip('\n')[9:151]
		R2_seq_plus=linecache.getline(R2,i-1).strip('\n')
		R2_MB_qua=linecache.getline(R2,i).strip('\n')[0:8]
		R2_true_seq_qua=linecache.getline(R2,i).strip('\n')[9:151]
		
		for keys in MB.keys():
			for item in MB[keys]:
				if re.match(item,R1_MB):
					
					R1_MB=str(keys)
					
					for keys in MB.keys():
						for item in MB[keys]:
							if re.match(item,R2_MB):
								
								R2_MB=str(keys)
								
			
		if R1_MB in MB.keys() and R2_MB in MB.keys():
			R1_seq_ID=str(R1_seq_ID+':'+str(R1_MB)+'+'+str(R2_MB))
			R2_seq_ID=str(R2_seq_ID+':'+str(R1_MB)+'+'+str(R2_MB))
			R1_output.write(R1_seq_ID+'\n'+R1_true_seq+'\n'+R1_seq_plus+'\n'+R1_true_seq_qua+'\n')
			R2_output.write(R2_seq_ID+'\n'+R2_true_seq+'\n'+R2_seq_plus+'\n'+R2_true_seq_qua+'\n')
	
		
	
	i+=1
	if i%4000000==0:
		print (str(int(i/4000000))+'M...reads processed')	
R1_output.close()
R2_output.close()

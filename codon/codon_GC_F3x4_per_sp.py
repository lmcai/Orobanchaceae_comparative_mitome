from Bio import SeqIO
import os
files=os.listdir('./')
files=[i for i in files if i.endswith('.codon.fas')]

for file in files:
	recs=SeqIO.parse(file,'fasta')
	out=open(file.split('.')[0]+'.GCperGene.tsv','a')
	for rec in recs:
		seq_str=str(rec.seq)
		proceed=len(seq_str) % 3
		if proceed==0:
			gc1=0
			gc2=0
			gc3=0
			for j in range(0,len(seq_str)):
				codon=j%3
				if codon==0:
					if seq_str[j]=='C' or seq_str[j]=='G' or seq_str[j]=='c' or seq_str[j]=='g':gc1=gc1+1
				elif codon==1:
					if seq_str[j]=='C' or seq_str[j]=='G' or seq_str[j]=='c' or seq_str[j]=='g':gc2=gc2+1
				elif codon==2:
					if seq_str[j]=='C' or seq_str[j]=='G' or seq_str[j]=='c' or seq_str[j]=='g':gc3=gc3+1
			d=out.write('\t'.join([rec.id,str(gc1),str(gc2),str(gc3),str(len(seq_str))])+'\n')
		else:
			print(file,rec.id)
	out.close()
	
					

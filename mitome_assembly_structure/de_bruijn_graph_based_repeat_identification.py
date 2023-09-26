import os
from numpy import median

#calculate average kmer coverage
files = os.listdir('./')
files=[i for i in files if i.endswith('gfa')]
out=open('average_kmer_cov.tsv','a')
for file in files:
	x=open(file).readlines()
	a=[]
	for l in x:
		if l.startswith('S'):
			LN=l.split()[3]
			LN=int(LN.split(':')[-1])
			LC=l.split()[4]
			LC=float(LC.split(':')[-1])
			a.append(LC/LN)
	out.write(file+'\t'+str(median(a))+'\n')

out.close()

#identify repeat in raw graph
import os

files = os.listdir('./')
files=[i for i in files if i.endswith('gfa')]
out=open('repeat_candidate.tsv','a')
for file in files:
	x=open(file).readlines()
	contig_len={}
	contig_cov={}
	contig_connectness={}
	all_cov=[]
	for l in x:
		if l.startswith('S'):
			LN=l.split()[3]
			LN=int(LN.split(':')[-1])
			LC=l.split()[4]
			LC=float(LC.split(':')[-1])
			contig_len[l.split()[1]]=LN
			contig_cov=LC/LN
			all_cov.append(LC/LN)
		elif l.startswith('L'):
			try:contig_connectness[l.split()[1]]=contig_connectness[l.split()[1]]+1
			except KeyError:contig_connectness[l.split()[1]]=1
			try:contig_connectness[l.split()[3]]=contig_connectness[l.split()[3]]+1
			except KeyError:contig_connectness[l.split()[3]]=1
	for k in contig_connectness.keys:
		if contig_connectness[k]>3:
			out.write(file+'\t'+k+'\t'+contig_len[k]+'\t'+str(contig_cov[k])+'\n')

out.close()
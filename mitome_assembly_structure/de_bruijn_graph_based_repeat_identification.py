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
from nympy import mean

files = os.listdir('./')
files=[i for i in files if i.endswith('gfa')]
out=open('repeat_candidate.tsv','a')
out.write('Species\tcontig_ID\tlength\tcoverage\tmean_neighbour_coverage\n')
for file in files:
	x=open(file).readlines()
	y=open(file.split('.gfa')[0]+'.scaffold.list').readlines()
	filtered_contigs=[i.strip() for i in y]
	filtered_contigs=list(set(filtered_contigs))
	contig_len={}
	contig_cov={}
	contig_connectness={}
	contig_neighbours={}
	#contig_seq={}
	for l in x:
		if l.startswith('S'):
			#contig_seq[l.split()[1]]=l.split()[2]
			LN=l.split()[3]
			LN=int(LN.split(':')[-1])
			LC=l.split()[4]
			LC=float(LC.split(':')[-1])
			contig_len[l.split()[1]]=LN
			contig_cov[l.split()[1]]=LC/LN
		elif l.startswith('L'):
			try:
				contig_connectness[l.split()[1]]=contig_connectness[l.split()[1]]+1
				contig_neighbours[l.split()[1]].append(l.split()[3])
			except KeyError:
				contig_connectness[l.split()[1]]=1
				contig_neighbours[l.split()[1]]=[l.split()[3]]
			try:
				contig_connectness[l.split()[3]]=contig_connectness[l.split()[3]]+1
				contig_neighbours[l.split()[3]].append(l.split()[1])
			except KeyError:
				contig_connectness[l.split()[3]]=1
				contig_neighbours[l.split()[3]]=[l.split()[1]]
	#output coverage~length to tsv
	out1=open(file.split('.gfa')[0]+'.cov_len.tsv','a')
	for k in contig_connectness.keys():
		#write coverage~length spreadsheet for all contigs
		d=out1.write(k+'\t'+str(contig_cov[k])+'\t'+str(contig_len[k])+'\n')
		#check if this contig is in the final assembly
		in_final_assembly=0
		contig_elements=k.split('_')
		if len([j for j in contig_elements if j in filtered_contigs])>0:in_final_assembly=1
		#calculate mean kmer coverager for all neighbors
		mean_neighbour_cov=[str(int(contig_cov[j])) for j in contig_neighbours[k]]
		if contig_connectness[k]>3 and in_final_assembly==1:
			out.write(file+'\t'+k+'\t'+str(contig_len[k])+'\t'+str(int(contig_cov[k]))+'\t'+str('_'.join(mean_neighbour_cov))+'\n')
	out1.close()


out.close()
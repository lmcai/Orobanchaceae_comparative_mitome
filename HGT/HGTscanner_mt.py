from Bio import SeqIO
import os,argparse
import datetime
from ete3 import Tree
import ete3
import pybedtools
from numpy import median
parser = argparse.ArgumentParser(description='HGTscanner_mt is a utility wrappepr to identify HGT blocks in organellar (mostly mitochondrial) genomes.')
parser.add_argument('-q', metavar='query', help='fasta file of the target mitome', required=True)
parser.add_argument('-ref', metavar='reference', help='one fasta file containing all custom references of both close relatives and potential HGT donor. This will be combined with the NCBI Viridiplantae mito database.')
parser.add_argument('-o', metavar='output', help='output prefix', required=True)
parser.add_argument('-f', metavar='family', help='family of the query for HGT classification', required=True)
parser.add_argument('-b', metavar='mask', help='bed file for regions to be masked')

args = parser.parse_args()


query=args.q
if args.ref:
	reference=args.ref
sp=args.o
fam=args.family

###############
#MASK gene/MTPT regions
if args.b:
	print(str(datetime.datetime.now())+'\tMasking query sequence '+query+' using the bed file: '+args.b)

###########
#BLAST
if args.ref:
	#add custom references to database
	print(str(datetime.datetime.now())+'\tAdded custom reference '+reference+' to the NCBI Viridiplantae mitochondrial database')
	S='cat '+reference+' Viridiplantae_mt.fasta >mt_db.fas'
	os.system(S)
	S='makeblastdb -in '+query+' -out mt -dbtype nucl >/dev/null'
	os.system(S)
	S='blastn -task dc-megablast -query '+query+' -db mt -outfmt 6 -evalue 1e-20 >'+sp+'.temp.blast'
	os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
else:
	#use built-in database
	S='blastn -task dc-megablast -query '+query+' -db mt -outfmt 6 -evalue 1e-20 >'+sp+'.temp.blast'
	os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)


####################################
#Add taxonomic information to each hit at species and family level
x=open(sp+'.temp.blast').readlines()
out=open(sp+'.taxon.blast','a')
species={}
family={}
y=open('Viridiplantae_mt.taxonomy').readlines()
for l in y:
	species[l.split('\t')[0]]=l.split('\t')[1]
	family[l.split('\t')[0]]=l.split('\t')[2].strip()


for l in x:
	try:
		d=out.write(l.strip()+'\t'+species[l.split()[1]]+'\t'+family[l.split()[1]]+'\n')
	except KeyError:
		d=out.write(l.strip()+'\t\t\n')


out.close()

###########
#sort blast results, give each row a unique id
S="cat "+sp+".taxon.blast | sort -k1,1 -k7,7n -k11,11n | awk -v OFS='\\t' '{if ($8-$7 <= 20000) print $1, $7, $8, $2, $9, $10, $13, $14, NR}' > "+sp+".taxon_sorted.bed"
os.system(S)

####################################
#define potential HGT blocks
S="bedtools merge -i "+sp+".taxon_sorted.bed -c 9 -o collapse >"+sp+'.temp.bed'
os.system(S)
#blocks <1500 bp in length can be used directly for downstream HGT inference; larger blocks >1500bp, break it into smaller ones because it can be caused by blasting to close relative or even itself
x=open(sp+'.temp.bed').readlines()
y=open(sp+'.taxon_sorted.bed').readlines()
out=open(sp+'.hgt.bed','a')

def id2bed(ids,bed_file):
	id_dict={}
	for l in bed_file:
		id_dict[l.split()[-1]]=l
	filtered_bed=[id_dict[j] for j in ids]
	return(filtered_bed)

#**the most important function in this pipeline**
def break_large_beds(bed_file):
	#remove blast hits from the same family to build new merged bed file (hopefully they are smaller). These hits will be added back don't worry
	filtered=[]
	same_fam=[]
	for j in bed_file:
		#replace empty columns with 'NA', pybedtools complain about them
		if j.split('\t')[6]=='':
			j='\t'.join(j.split('\t')[:6])+'\t'+j.split()[3]+'\t'+'\t'.join(j.split('\t')[7:])
		if j.split('\t')[7]=='':
			j='\t'.join(j.split('\t')[:7])+'\tNA\t'+'\t'.join(j.split('\t')[8:])
		if not (j.split('\t')[7]==fam or j.split('\t')[3].startswith('Oro')):
			filtered.append(j)
		else:
			same_fam.append(j)
	filtered_bed = pybedtools.BedTool(''.join(filtered), from_string=True)
	samefam_bed = pybedtools.BedTool(''.join(same_fam), from_string=True)
	merged_bed = filtered_bed.merge(c=9,o='collapse')
	#if merged_bed meet the criteria, you can define interactions with '+' and '-'!!
	if ###:
		(merged_bed+samefam_bed).count()
	else:
		#iteratively remove the longest until meet the criteria?


for l in x:
	ids=l.split('\t')[-1].strip()
	raw_beds=id2bed(ids.split(','),y)
	median_len=median([int(j.split()[2])-int(j.split()[1]) for j in raw_beds])
	bed_len=int(l.split()[2])-int(l.split()[1])
	#get median blast hit length, if <500 bp or largely overlap with the bed size, output
	if median_len>0.5*bed_len or bed_len<500:
		out.write(l)
	else:
		new_beds=break_large_beds(raw_beds)
		
		
		
		















print(str(datetime.datetime.now())+'\tFound '+str()+' homologous genetic blocks for further exaination')

####################################
#extract sequences for each block
loci=open(sp+'.temp.bed').readlines()
hits=open(sp+'.mtpt.bed').readlines()
q_recs=SeqIO.index(query,'fasta')
if args.ref:
	ref_recs=SeqIO.index('mt_db.fas', 'fasta')
else:
	ref_recs=SeqIO.index('Viridiplantae_mt.fasta', 'fasta')

seq_loc={}
for l in hits:
	seq_loc[l.split()[-1]]=l

for l in loci:
	out=open(sp+'.temp.'+str(order)+'.fas','w')
	d=SeqIO.write(q_recs[l.split()[0]][(int(l.split()[1])-1):int(l.split()[2])],out,'fasta')
	out.close()
	out=open(sp+'.temp.'+str(order)+'.fas','a')
	bed_ids=l.split()[3].strip()
	bed_ids=bed_ids.split(',')
	for i in bed_ids:
		line=seq_loc[i]
		id=line.split()[3]
		start=int(line.split()[4])-1
		end=int(line.split()[5])
		d=out.write('>'+id+'_'+str(start)+'_'+str(end)+'\n')
		d=out.write(str(ref_recs[id].seq[start:end])+'\n')
	out.close()
	order=order+1

print(str(datetime.datetime.now())+'\tExatracted potential HGT sequences for '+sp)

#############
#alignment and phylogenetic reconstruction
print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')

for i in range(1,order):
	S="mafft --quiet --adjustdirection "+sp+".temp."+str(i)+".fas | sed 's/_R_//g' > "+sp+".gt."+str(i)+".aln.fas"
	os.system(S)
	b=SeqIO.index(sp+".gt."+str(i)+".aln.fas",'fasta')
	q=loci[i-1].split()[0]
	if len(b[q].seq)<10000:
		S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".gt."+str(i)+".aln.fas >/dev/null 2>&1"
		os.system(S)
		print(str(datetime.datetime.now())+'\tLoci #'+str(i))
	else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')

os.system('rm '+sp+'*.bionj')
os.system('rm '+sp+'*.gz')
os.system('rm '+sp+'*.log')
os.system('rm '+sp+'*.iqtree')
os.system('rm '+sp+'*.mldist')
os.system('rm '+sp+'*.phy')
os.system('rm '+sp+'*.contree')
os.system('rm '+sp+'*.nex')

##############
#Evaluate the source of the region and output summary file
out=open(sp+'.hgt.sum.tsv','w')
out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
out.close()
out=open(sp+'.hgt.sum.tsv','a')

for i in range(1,order):
	q=loci[i-1].split()[0]
	outgroup=[]
	try:
		t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
		ancestor=t.get_midpoint_outgroup()
		t.set_outgroup(ancestor)
		q_branch=t&q
		if not q_branch.get_ancestors()[0].is_root():
			sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
		else:
			t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
			outgroup=[leaf.name for leaf in t if leaf.name.startswith('Sorghum')]
			if len(outgroup)>0:
				t.set_outgroup(outgroup[0])
				q_branch=t&q
				sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
			else:
				outgroup=[leaf.name for leaf in t if leaf.name.startswith('Rehm')]
				if outgroup:
					t.set_outgroup(t&outgroup[0])
					q_branch=t&q
					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
				else:
					#no sorghum no rehmannia
					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
	except ete3.parser.newick.NewickError:
		sisters=open(sp+'.temp.'+str(i)+".fas").readlines()
		sisters=[i[1:].strip() for i in sisters if (i.startswith('>')) and (not i[1:].strip()==q)]
		out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+'ALL HITS: '+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
		
out.close()
os.system('rm '+sp+'.temp.*.fas')
os.system('rm '+sp+'.temp.bed')
os.system('rm '+sp+'.n*')

if not os.path.isdir('HGTscanner_supporting_files'):os.mkdir('HGTscanner_supporting_files')
os.system('mv '+sp+'.*.aln.fas* HGTscanner_supporting_files')
print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See summary file in '+sp+'.hgt.sum.tsv')
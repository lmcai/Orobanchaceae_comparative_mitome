from Bio import SeqIO
import os,argparse
import datetime
from ete3 import Tree
import ete3
import pybedtools
from numpy import median
from statistics import mode

parser = argparse.ArgumentParser(description='HGTscanner_mt is a utility wrappepr to identify HGT blocks in organellar (mostly mitochondrial) genomes.')
parser.add_argument('-q', metavar='query', help='fasta file of the target mitome', required=True)
parser.add_argument('-ref', metavar='reference', help='one fasta file containing all custom references of both close relatives and potential HGT donor. This will be combined with the NCBI Viridiplantae mito database.')
parser.add_argument('-o', metavar='output', help='output prefix', required=True)
parser.add_argument('-f', metavar='family', help='family of the query for HGT classification', required=True)
parser.add_argument('-b', metavar='mask', help='bed file for regions to be masked')

#pass argument values to variables
args = parser.parse_args()
query=args.q
if args.ref:
	reference=args.ref
sp=args.o
fam=args.f
if args.b:
	mask_bed=args.b


###############
#MASK gene/MTPT regions
def mask_fasta_with_bed(fasta_file, bed_file, output_file):
	sequences = SeqIO.index(fasta_file, "fasta")
	bed = open(bed_file).readlines()
	recs={}
	for k in sequences.keys():
		recs[k]=str(sequences[k].seq)
	for l in bed:
   		sequence_id = l.split()[0]
   		start = int(l.split()[1])
   		end = int(l.split()[2])
   		sequence = recs[sequence_id]
   		masked_sequence = sequence[:(start-1)] + "N" * (end - start + 1) + sequence[end:]
   		recs[sequence_id] = masked_sequence
	out=open(output_file,'w')
	for k in recs.keys():
		d=out.write('>'+k+'\n'+recs[k]+'\n')



if args.b:
	print(str(datetime.datetime.now())+'\tMasking query sequence '+query+' using the bed file: '+args.b)
	mask_fasta_with_bed(query, mask_bed, sp+'.mask.fas')


###########
#BLAST
if args.ref:
	#add custom references to database
	print(str(datetime.datetime.now())+'\tAdded custom reference '+reference+' to the NCBI Viridiplantae mitochondrial database')
	S='cat '+reference+' Viridiplantae_mt.fasta >'+sp+'.mt_db.fas'
	os.system(S)
else:
	S='cp Viridiplantae_mt.fasta '+sp+'.mt_db.fas'
	os.system(S)


S='makeblastdb -in '+sp+'.mt_db.fas -out '+sp+'.mt -dbtype nucl >/dev/null'
os.system(S)


if args.b:S='blastn -task dc-megablast -query '+sp+'.mask.fas -db '+sp+'.mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
else:S='blastn -task dc-megablast -query '+query+' -db '+sp+'.mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
os.system(S)
print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)


####################################
#Add taxonomic information to each hit at species and family level
x=open(sp+'.raw.blast').readlines()
out=open(sp+'.taxon.blast','w')
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
		if l.split()[1].startswith('Oro'):d=out.write(l.strip()+'\t'+l.split()[1]+'\tOrobanchaceae\n')
		else:d=out.write(l.strip()+'\t'+l.split()[1]+'\tNA\n')


out.close()

###########
#sort blast results, give each row a unique id, only hits <20k are printed
S="cat "+sp+".taxon.blast | sort -k1,1 -k7,7n -k11,11n | awk -v OFS='\\t' '{if ($8-$7 <= 20000) print $1, $7, $8, $2, $9, $10, $13, $14, NR}' > "+sp+".taxon_sorted.bed"
os.system(S)

####################################
#define potential HGT blocks

#classify these hits based on source, here, alignments from Orobanchaceae, Lamiaceae, Oleaceae, and Scrophulariaceae will be regarded as close relatives
#This is necessary because they create long synteny blocks that may be broken by HGT in the middle

x=open(sp+".taxon_sorted.bed").readlines()
otherfam=[]
samefam=[]
for l in x:
	if not l.split()[7] in ['Orobanchaceae','Lamiaceae','Scrophulariaceae','Oleaceae','Phrymaceae','Acanthaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Bignoniaceae','Anacardiaceae']:
		otherfam.append(l)
	else:
		samefam.append(l)

otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=9,o='collapse')
samefam_bed=pybedtools.BedTool(''.join(samefam), from_string=True)
out=open(sp+'.merged.bed','w')
d=out.write(str(otherfam_merged))

print(str(datetime.datetime.now())+'\tFound '+str(len(otherfam_merged))+' homologous genetic blocks for further examination')

###############everything below in this section is abandoned
#def bedOK2go(bed_txt,all_aln):
#	ids=bed_txt.split('\t')[-1].strip()
#	aln=id2bed(ids.split(','),all_aln)
#	median_len=median([int(j.split()[2])-int(j.split()[1]) for j in aln])
#	bed_start=min([int(j.split()[1]) for j in aln])
#	bed_end=max([int(j.split()[2]) for j in aln])
#	bed_len=bed_end - bed_start
	#get median blast hit length, if <500 bp or largely overlap with the bed size, output
#	if median_len>0.5*bed_len or bed_len<500:return(1)
#	else:return(0)

#def break_up_range(start, end, step):
#    ranges = []
#    for i in range(start, end + 1, step):
#        ranges.append((i, min(i + step - 1, end)))
#    return ranges


#This function is abandoned since I did not find a good way to cluster individual hits to groups with minor differences on their boundries
#def break_large_beds(bed_file):
	#I. remove blast hits from the same family to build new merged bed file (hopefully they are smaller). These hits will be added back don't worry
#	filtered=[]
#	addback=[]
#	for j in bed_file:
#		if not (j.split('\t')[7]==fam or j.split('\t')[3].startswith('Oro')):
#			filtered.append(j)
#		else:
#			addback.append(j)
#	filtered_bed = pybedtools.BedTool(''.join(filtered), from_string=True).merge(c=9,o='collapse')
	#if merged_bed meet the criteria, you can define interactions with '+' and '-'!!
#	filtered_bed_str=str(filtered_bed).strip()
#	OK=1
#	for l in filtered_bed_str.split('\n'):
#		if not bedOK2go(l,y):OK=0
#	if OK:
		#removing same family aln is sufficient to create acceptable homologs
#		addback_ids=[j.split()[-1] for j in addback]
		#add back aln ids to each smaller bed region 
#		filtered_bed_str=filtered_bed_str.strip()
#		output_txt=''
#		for l in filtered_bed_str.split('\n'):
#			output_txt=output_txt+l+','+','.join(addback_ids)+'\n'
#		return(output_txt)
#	else:
		#II. break up blocks based on coverage patterns
#		return('')
		#full_bed=str(pybedtools.BedTool(''.join(bed_file), from_string=True).merge())
		#broken_range=break_up_range(int(full_bed.split()[1]),int(full_bed.split()[2]),5)
		#broken_bed_txt=[full_bed.split()[0]+'\t'+str(i[0])+'\t'+str(i[1]) for i in broken_range]
		#broken_bed=pybedtools.BedTool('\n'.join(broken_bed_txt), from_string=True)
		#all_bed=pybedtools.BedTool(''.join(bed_file), from_string=True)
		#bed_cov=broken_bed.coverage(all_bed)
		
		#filtered_bed_len={}
		#filtered_bed_len_list=[]
		#for j in filtered:
		#	try:
		#		filtered_bed_len[str(int(j.split()[2])-int(j.split()[1]))].append(j)
		#	except KeyError:filtered_bed_len[str(int(j.split()[2])-int(j.split()[1]))]=[j]
		#	filtered_bed_len_list.append(int(j.split()[2])-int(j.split()[1]))
		#while True:
		#	cur_max=max(filtered_bed_len_list)
		#	cur_max_bed=filtered_bed_len[cur_max]
		#	for k in cur_max_bed:
		#		filtered.remove(k)
		#		addback.append(k)
		#		filtered_bed_len_list.remove(cur_max)
		#		filtered_bed = pybedtools.BedTool(''.join(filtered), from_string=True).merge()
		#		if bedOK2go(str(filtered_bed)):
					#save
		#			break

####################################
#extract sequences for each block

q_recs=SeqIO.index(query,'fasta')
ref_recs=SeqIO.index(sp+'.mt_db.fas', 'fasta')

def id2bed(ids,bed_file):
	id_dict={}
	for l in bed_file:
		id_dict[l.split()[-1]]=l
	filtered_bed=[id_dict[j] for j in ids]
	return(filtered_bed)

#takes in bed file of aligned region, order them, and filter for ones selected to output
#e.g.
#Pch_1	133169	133918	BA000042.1	329847	330607	Nicotiana_tabacum
#Pch_1	133907	134389	BA000042.1	330960	331502	Nicotiana_tabacum
#Pch_1	134385	136416	BA000042.1	331627	333652	Nicotiana_tabacum
#Pch_1	133169	133912	CP129458.1	165029	164277	Quercus_variabilis
#Pch_1	133169	133299	CP129458.1	220967	221096	Quercus_variabilis #this row will be removed
#Pch_1	133907	135460	CP129458.1	163913	162355	Quercus_variabilis
def aln_scaffolder(bedtxt):
	rawtable = [line.split('\t') for line in bedtxt]
	sortedtable = sorted(rawtable, key=lambda x: (x[3], int(x[1])))
	filteredtable=[]
	cursp=''
	end_pos=0
	for l in sortedtable:
		if l[3]!=cursp:
			#a new sp
			filteredtable.append(l)
			cursp=l[3]
			end_pos=int(l[2])
		else:#same sp
			if int(l[1])>end_pos-50:
				#allow for 50 bp overlap in the reference
				filteredtable.append(l)
				end_pos=int(l[2])
	consolidatedtable={}
	consolidated_block_count={}
	for l in filteredtable:
		try:
			consolidatedtable[l[3]].append(l)
			consolidated_block_count[l[3]]=consolidated_block_count[l[3]]+1
		except KeyError:
			consolidatedtable[l[3]]=[l]
			consolidated_block_count[l[3]]=1
	mode_consolidated_block_count=mode([consolidated_block_count[k] for k in consolidated_block_count.keys()])
	if mode_consolidated_block_count<2:
		outputtable=[]
		for k in consolidatedtable.keys():
			refpos=[l[1]+'-'+l[2] for l in consolidatedtable[k]]
			targetpos=[l[4]+'-'+l[5] for l in consolidatedtable[k]]
			outputtable.append(consolidatedtable[k][0][0]+'\t'+';'.join(refpos)+'\t'+consolidatedtable[k][0][3]+'\t'+';'.join(targetpos)+'\t'+'\t'.join(consolidatedtable[k][0][6:]))
	else:
		#break this scaffold into smaller chunks
		for k in consolidated_block_count.keys():
			outputtable=['X']
			if consolidated_block_count[k]==mode_consolidated_block_count:
				outputtable.append(['\t'.join(i[:3]) for i in consolidatedtable[k]])
				break
	return(outputtable)

def seq2seq_ortho_extraction(seed_file,targets_file,output_handle):
	S='makeblastdb -in '+targets_file+' -out '+sp+'.temp -dbtype nucl >/dev/null'
	os.system(S)
	S='blastn -task dc-megablast -query '+seed_file+' -db '+sp+'.temp -outfmt 6 -evalue 1e-20 | sort -k2,2 -k11,11n> '+sp+'.temp.blast'
	os.system(S)
	hits=open(sp+'.temp.blast').readlines()
	#select only the best hit per target file
	cur_hit=''
	for l in hits:
		ncbiID=l.split()[1]
		if ncbiID!=cur_hit:
			try:
				d=out.write('>'+family[ncbiID]+'|'+ncbiID+'|'+species[ncbiID]+'\n')
			except KeyError:
				if ncbiID.startswith('Oro'):d=out.write('>Orobanchaceae|'+ncbiID+'\n')
				else:d=out.write('>NA|'+ncbiID+'\n')
			if int(l.split()[8])<int(l.split()[9]):
				d=out.write(str(ref_recs[l.split()[1]].seq[(int(l.split()[8])-1):int(l.split()[9])])+'\n')
			else:
				d=out.write(str(ref_recs[l.split()[1]].seq[(int(l.split()[9])-1):int(l.split()[8])])+'\n')
			cur_hit=ncbiID
		else:
			pass

order=1
num=1
mapout=open(sp+'.alnmap.bed','w')
for hit in otherfam_merged:
	#gather overlapping alignment from both other families and close relatives
	ids=hit.fields[3]
	raw_beds_txt=id2bed(ids.split(','),x)
	seqout_beds=aln_scaffolder(raw_beds_txt)
	if seqout_beds[0]=='X':
		#this hit needs to be further divided into smaller chunks
		raw_beds=pybedtools.BedTool(''.join(raw_beds_txt), from_string=True)
		for i in seqout_beds[1]:
			out=open(sp+'.hgt.'+str(order)+'.fas','w')
			subhit=pybedtools.BedTool(i, from_string=True)
			new_raw_beds=raw_beds.intersect(subhit,wa=True,f=0.4)
			#write other family
			min_start=1000000
			max_end=0
			for l in new_raw_beds:
				l=str(l).split()
				d=out.write('>'+l[7]+'|'+l[3]+'|'+l[4]+'-'+l[5]+'|'+l[6]+'\n')
				start=int(l[4])
				end=int(l[5])
				if int(l[1])<min_start:min_start=int(l[1])
				if int(l[2])>max_end:max_end=int(l[2])
				if start<end:
					d=out.write(str(ref_recs[l[3]].seq[(start-1):end])+'\n')
				else:
					d=out.write(str(ref_recs[l[3]].seq[(end-1):start].reverse_complement())+'\n')
			#write query
			d=SeqIO.write(q_recs[l[0]][(min_start-1):max_end],out,'fasta')
			#write close relative
			samefam_hit=samefam_bed.intersect(subhit,f=0.4)
			if str(samefam_hit)!='':
				d=SeqIO.write(q_recs[l[0]][(min_start-1):max_end],sp+'.tempseed.fas','fasta')
				out2=open(sp+'.tempTarget.fas','w')
				for ll in samefam_hit:
					d=SeqIO.write(ref_recs[ll.fields[3]],out2,'fasta')
				out2.close()
				seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas',out)
			d=mapout.write(hit.chrom+'\t'+str(min_start)+'\t'+str(max_end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
			order=order+1
	else:
		out=open(sp+'.hgt.'+str(order)+'.fas','w')
		#write other family
		for l in seqout_beds:
			d=out.write('>'+'|'.join([l.split()[5]]+l.split()[2:5])+'\n')
			secs=l.split()[3]
			sequence=''
			for sec in secs.split(';'):
				start=int(sec.split('-')[0])
				end=int(sec.split('-')[1])
				if start<end:
					sequence=sequence+str(ref_recs[l.split()[2]].seq[(start-1):end])	
				else:
					sequence=sequence+str(ref_recs[l.split()[2]].seq[(end-1):start].reverse_complement())
				d=out.write(sequence+'\n')	
		#write query
		d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],out,'fasta')
		#write close relatives
		hit_bed=pybedtools.BedTool(str(hit),from_string=True)
		samefam_hit=samefam_bed.intersect(hit_bed)
		d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],sp+'.tempseed.fas','fasta')
		out2=open(sp+'.tempTarget.fas','w')
		for l in samefam_hit:
			d=SeqIO.write(ref_recs[l.fields[3]],out2,'fasta')
		out2.close()
		seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas',out)
		d=mapout.write(hit.chrom+'\t'+str(hit.start)+'\t'+str(hit.end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
		order=order+1
		out.close()
	current_time = datetime.datetime.now()
	print(f"{current_time}\tExtracting sequences from homologous genetic block #{num}", end='\r')
	num=num+1
		
mapout.close()		
print(f"{current_time}\tA total of #{order} aligned sequences from #{num} merged homologous genetic blocks were extracted.", end='\r')

#############
#alignment and phylogenetic reconstruction
print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')

#for i in range(1,order):
#	current_time = datetime.datetime.now()
#	print(f"{current_time}\t Sequence alignment and IQTREE for alignment #{i}", end='\r')
#	S="timeout 20m mafft --genafpair --maxiterate 1000 --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
#	os.system(S)
	#recs=list(SeqIO.parse(sp+".hgt."+str(i)+".fas",'fasta'))
	#max_len=max([len(rec.seq) for rec in recs])
	#if max_len<1500:
	#	S="mafft --localpair --maxiterate 1000 --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
	#	S="mafft --adjustdirection --6merpair --addfragments othersequences referencesequence > output"
	#	os.system(S)
	#	S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".hgt."+str(i)+".aln.fas >/dev/null 2>&1"
	#	os.system(S)
	#elif max_len<5000:
	#	S="mafft --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
	#	os.system(S)
	#	S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".hgt."+str(i)+".aln.fas >/dev/null 2>&1"
	#	os.system(S)
	#else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')

#os.system('rm '+sp+'*.bionj')
#os.system('rm '+sp+'*.gz')
#os.system('rm '+sp+'*.log')
#os.system('rm '+sp+'*.iqtree')
#os.system('rm '+sp+'*.mldist')
#os.system('rm '+sp+'*.phy')
#os.system('rm '+sp+'*.contree')
#os.system('rm '+sp+'*.nex')

os.system('rm '+sp+'.raw.blast')
os.system('rm '+sp+'.temp*')
os.system('rm '+sp+'.mt.n*')
os.system('rm '+sp+'.mt_db.fas')

if not os.path.isdir(sp+'_HGTscanner_supporting_files'):os.mkdir(sp+'_HGTscanner_supporting_files')
os.system('mv '+sp+'.hgt.*.fas '+sp+'_HGTscanner_supporting_files')
print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See sequence file in '+sp+'_HGTscanner_supporting_files')
#parse codeml output for codon counts and GC1, GC2, GC3
#first run an codeml analysis with one omega for the entire tree
import re
x=open('sp43_codeml_1omega.mlc').readlines()
codon_counts={}
#locate codon count block
codon=[]
for l in x[6:]:
	#end of the block
	if l.startswith('Codon'):break
	else:
		if not (l.startswith('-') or l.strip() == ""):
			temp=l.split('|')
			for i in temp:
				matches = re.match(r'([a-zA-Z\s\*]+)([\d\s]+)', i)
				try:codon_counts[matches.group(1).strip()]=codon_counts[matches.group(1).strip()]+matches.group(2).strip().split()
				except KeyError:
					codon_counts[matches.group(1).strip()]=matches.group(2).strip().split()
					codon.append(matches.group(1).strip())
	
#parse the codon gc count block
sp_all=[]
sp='Species'
gc1='GC1'
gc2='GC2'
gc3='GC3'
out=open('GC123.tsv','a')
for l in x[10:]:
	if l.startswith('#'):
		#output to file
		d=out.write(sp+'\t'+str(gc1)+'\t'+str(gc2)+'\t'+str(gc3)+'\n')
		#start the new round
		sp=l.split()[1]
		sp_all.append(sp)
	elif l.startswith('position  1'):
		temp=l.split()
		gc1=float(temp[3].split(':')[1])+float(temp[5].split(':')[1])
	elif l.startswith('position  2'):
		temp=l.split()
		gc2=float(temp[3].split(':')[1])+float(temp[5].split(':')[1])
	elif l.startswith('position  3'):
		temp=l.split()
		gc3=float(temp[3].split(':')[1])+float(temp[5].split(':')[1])
	#start of gc block
	elif l.startswith('Sum'):break

#output the final species

d=out.write(sp+'\t'+str(gc1)+'\t'+str(gc2)+'\t'+str(gc3)+'\n')
out.close()
#assuming the codon and the gc data are in the same order, output codon table with species names
out=open('codon_table.tsv','a')
out.write('Codon\t'+'\t'.join(sp_all)+'\n')
for i in codon:
	d=out.write(i+'\t'+'\t'.join(codon_counts[i])+'\n')

out.close()

	
import os
from Bio import SeqIO

files=os.listdir('./')
out=open('IR_len_sum.tsv','a')
for file in files:
	REPs=[]
	recs=SeqIO.read(file,'genbank')
	total=int(recs.features[0].location.end)
	# Initialize variables to store the positions of the two gene features
	gene1_start = None
	gene1_end = None
	gene2_start = None
	gene2_end = None
	for feature in recs.features:
		if feature.type == "repeat_region":
			REPs.append(feature)
			
	if len(REPs)==2:
		gene1_start = REPs[0].location.start
		gene1_end = REPs[0].location.end
		gene2_start = REPs[1].location.start
		gene2_end = REPs[1].location.end

	# Check if both gene features were found
	if gene1_start is not None and gene2_start is not None:
    	# Calculate the distance between the two gene features
    	SR1 = abs(int(gene2_start) - int(gene1_end))
    	IR = abs(int(gene1_end) - int(gene1_start))
		SR2 = total-SR1-2*IR
		SSR=min([SR1,SR2])
		LSR=max([SR1,SR2])
		d=out.write(file+'\t'+str(IR)+'\t'+str(SSR)+'\t'+str(LSR)+'\n')
	else:
		d=out.write(file+'\n')
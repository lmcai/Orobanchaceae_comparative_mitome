from Bio import SeqIO
#x=open('/Users/lcai/Dropbox/Parasitic_plant_mitonuclear_coevolution/Comparative_mitogenome/mt_gene.txt').readlines()
#x=[i.strip() for i in x]

#scaffolds
records = SeqIO.parse('/Users/lcai/Dropbox/Parasitic_plant_mitonuclear_coevolution/Organelle_assemblies/mt_assembly/Aphyllon_purpureum_sub3.gb', "genbank")
i=1
for rec in records:
	for feature in rec.features:
		if feature.type=='gene':
			gene_name=feature.qualifiers['gene'][0]
			outfile=open(gene_name+'.fas','a')
			d=outfile.write('>Monochasma_japonicum_sub'+str(i)+'\n')
			d=outfile.write(str(feature.extract(rec).seq)+'\n')
			outfile.close()
	i=i+1

#batch processing
files=os.listdir('./')
for f in files:
	try:
		records = SeqIO.parse(f, "genbank")
		i=1
		for rec in records:
			for feature in rec.features:
				if feature.type=='gene':
					gene_name=feature.qualifiers['gene'][0]
					outfile=open(gene_name+'.fas','a')
					d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'_sub'+str(i)+'\n')
					d=outfile.write(str(feature.extract(rec).seq)+'\n')
					outfile.close()
			i=i+1
	except:
		print(f)
		pass



#circularized
records = SeqIO.read('/Users/lcai/Dropbox/Parasitic_plant_mitonuclear_coevolution/Organelle_assemblies/mt_assembly/Aphyllon_purpureum_sub3.gb', "genbank")
for feature in records.features:
	if feature.type=='gene':
		gene_name=feature.qualifiers['gene'][0]
		outfile=open(gene_name+'.fas','a')
		d=outfile.write('>Orobanche_fasciculata\n')
		d=outfile.write(str(feature.extract(records).seq)+'\n')
		outfile.close()


##below is for fasta sequences
for rec in recs:
	g=rec.id
	out=open(g+'.fas','a')
	d=out.write('>Aeginetia_indica_NC069194\n')
	d=out.write(str(rec.seq)+'\n')
	out.close()
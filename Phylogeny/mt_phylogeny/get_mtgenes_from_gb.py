from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

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



#########################
#Get intronic regions of nad genes and rrn genes

files=os.listdir('./')
files=[f for f in files if f.endswith('gb')]

#get nad1,nad2,nad5 individual exons
for f in files:
	try:
		records = SeqIO.parse(f, "genbank")
		for rec in records:
			for feature in rec.features:
				if feature.type=='exon':
					try:gene_name=feature.qualifiers['gene'][0]
					except KeyError:pass
					exon_number=feature.qualifiers['number'][0]
					if (gene_name=='nad1' and exon_number=='1') or (gene_name=='nad1' and exon_number=='4') or (gene_name=='nad1' and exon_number=='5'):
						outfile=open(gene_name+'_exon'+exon_number+'.fas','a')
						f=f.split('.')[0]
						d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
						d=outfile.write(str(feature.extract(rec).seq)+'\n')
						outfile.close()
	except:
		print(f)
		pass

#get rrn18 and rrn26

for f in files:
	try:
		records = SeqIO.parse(f, "genbank")
		for rec in records:
			for feature in rec.features:
				#get rrn18 and rrn26
				if feature.type=='gene':
					gene_name=feature.qualifiers['gene'][0]
					if gene_name=='rrn18' or gene_name=='rrn26':
						outfile=open(gene_name+'.fas','a')
						f=f.split('.')[0]
						d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
						d=outfile.write(str(feature.extract(rec).seq)+'\n')
						outfile.close()
	except:
		print(f)
		pass

for f in files:
	try:
		records = SeqIO.parse(f, "genbank")
		#initiate locations of nad1 exon2-3; nad2 exon1-2, exon3-4-5; nad5 exon 1-2, exon 4-5
		for rec in records:
			nad1_exon23_locations = []
			nad2_exon12_locations = []
			nad2_exon345_locations = []
			nad5_exon12_locations = []
			nad5_exon45_locations = []
			for feature in rec.features:
				if feature.type=='exon':
					try:gene_name=feature.qualifiers['gene'][0]
					except KeyError:pass
					exon_number=feature.qualifiers['number'][0]
					#write multiple gene blocks, this involved combining regions using the 'sum' function
					if (gene_name=='nad1' and exon_number=='2') or (gene_name=='nad1' and exon_number=='3'):
						nad1_exon23_locations.append(feature.location)
					elif (gene_name=='nad2' and exon_number=='1') or (gene_name=='nad2' and exon_number=='2'):
						nad2_exon12_locations.append(feature.location)
					elif (gene_name=='nad2' and exon_number=='3') or (gene_name=='nad2' and exon_number=='5'):
						nad2_exon345_locations.append(feature.location)
					elif (gene_name=='nad5' and exon_number=='1') or (gene_name=='nad2' and exon_number=='2'):
						nad5_exon12_locations.append(feature.location)
					elif (gene_name=='nad5' and exon_number=='4') or (gene_name=='nad2' and exon_number=='5'):
						nad5_exon45_locations.append(feature.location)
			if len(nad1_exon23_locations)>1:
				sorted_nad1_exon23_locations = sorted(nad1_exon23_locations, key=lambda x: x.start)
				combined_nad1_exon23_location = FeatureLocation(sorted_nad1_exon23_locations[0].start, sorted_nad1_exon23_locations[-1].end)
				outfile=open('nad1_exon23.fas','a')
				f=f.split('.')[0]
				d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
				d=outfile.write(str(combined_nad1_exon23_location.extract(rec.seq))+'\n')
				outfile.close()
			if len(nad2_exon12_locations)>1:
				sorted_nad2_exon12_locations = sorted(nad2_exon12_locations, key=lambda x: x.start)
				combined_nad2_exon12_location = FeatureLocation(sorted_nad2_exon12_locations[0].start, sorted_nad2_exon12_locations[-1].end)
				outfile=open('nad2_exon12.fas','a')
				f=f.split('.')[0]
				d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
				d=outfile.write(str(combined_nad2_exon12_location.extract(rec.seq))+'\n')
				outfile.close()
			if len(nad2_exon345_locations)>1:
				sorted_nad2_exon345_locations = sorted(nad2_exon345_locations, key=lambda x: x.start)
				combined_nad2_exon345_location = FeatureLocation(sorted_nad2_exon345_locations[0].start, sorted_nad2_exon345_locations[-1].end)
				outfile=open('nad2_exon345.fas','a')
				f=f.split('.')[0]
				d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
				d=outfile.write(str(combined_nad2_exon345_location.extract(rec.seq))+'\n')
				outfile.close()
			if len(nad5_exon12_locations)>1:
				sorted_nad5_exon12_locations = sorted(nad5_exon12_locations, key=lambda x: x.start)
				combined_nad5_exon12_location = FeatureLocation(sorted_nad5_exon12_locations[0].start, sorted_nad5_exon12_locations[-1].end)
				outfile=open('nad5_exon12.fas','a')
				f=f.split('.')[0]
				d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
				d=outfile.write(str(combined_nad5_exon12_location.extract(rec.seq))+'\n')
				outfile.close()
			if len(nad5_exon45_locations)>1:
				sorted_nad5_exon45_locations = sorted(nad5_exon45_locations, key=lambda x: x.start)
				combined_nad5_exon45_location = FeatureLocation(sorted_nad5_exon45_locations[0].start, sorted_nad5_exon45_locations[-1].end)
				outfile=open('nad5_exon45.fas','a')
				f=f.split('.')[0]
				d=outfile.write('>'+f.split('_')[0]+'_'+f.split('_')[1]+'\n')
				d=outfile.write(str(combined_nad5_exon45_location.extract(rec.seq))+'\n')
				outfile.close()
	except:
		print(f)
		pass
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

def convert_genbank_to_gff(genbank_file, output_gff_file):
    records = SeqIO.parse(genbank_file, "genbank")
    with open(output_gff_file, "w") as output_gff:
        for record in records:
            for feature in record.features:
                if feature.type != "source":
                    gff_line = f"{record.name}\t{feature.type}\t{feature.qualifiers['gene'][0]}\t{feature.location.start + 1}\t{feature.location.end}\t.\t{feature.strand}\t.\tID={feature.qualifiers['gene'][0]}"
                    output_gff.write(gff_line + "\n")

files=os.listdir('./')
files=[f for f in files if f.endswith('.gb')]

for f in files:
	try:
		convert_genbank_to_gff(f,f.split('.gb')[0]+'.gff')
	except KeyError:print(f)
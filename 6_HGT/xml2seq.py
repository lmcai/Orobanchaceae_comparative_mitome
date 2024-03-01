from Bio.Blast import NCBIXML

def parse_blast_xml(xml_file):
	file_handle = open(xml_file)
	blast_records = NCBIXML.parse(file_handle)
	for blast_record in blast_records:#blast_record = all blast hits for one query sequence
		out=open(blast_record.query+'.fas','a')
		for alignment in blast_record.alignments:#multiple hits for one ncbi record
			header=alignment.accession+'_'+alignment.hit_def
			i=1
			for hsp in alignment.hsps:
				d=out.write('>'+header+'_'+str(i)+'\n'+hsp.sbjct+'\n')
				i=i+1
		out.close()

# Call the function to parse the BLAST XML and extract information
parse_blast_xml('M2GJXUTX016-Alignment.xml')

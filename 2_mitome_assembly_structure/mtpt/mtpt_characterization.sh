#blast identify plastid-like region
makeblastdb -in Castilleja_paramensis.plastome.fasta -dbtype nucl -out pt
blastn -db pt -query Castilleja_paramensis_mitome.fasta -outfmt 6 -evalue 1e-70 >Castilleja_paramensis.mtpt.blast

#Extract the following columns from blast output
#then use bedtools to merge blast hit and calculate total length
#1.  qseqid      query or source (gene) sequence id
#7.  qstart      start of alignment in query
#8.  qend        end of alignment in query
awk -v OFS='\t' '{print $1,$7,$8}' Castilleja_paramensis.mtpt.blast | sort -k1,1 -k2,2n >bedtools.in.bed
bedtools merge -i bedtools.in.bed >Castilleja_paramensis.mtpt.bed

#get sum
cat Castilleja_paramensis.mtpt.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'

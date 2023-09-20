makeblastdb -in mt_assembly/misc/$1.fasta -dbtype nucl -out mt
blastn -db mt -query pt_assembly/Rehmannia_glutinosa.fasta -outfmt 6 -evalue 1e-40 >$1.comparative_mtpt.blast
awk -v OFS='\t' '{print $1,$7,$8}' $1.comparative_mtpt.blast | sort -k1,1 -k2,2n >bedtools.in.bed
bedtools merge -i bedtools.in.bed >$1.comparative_mtpt.bed

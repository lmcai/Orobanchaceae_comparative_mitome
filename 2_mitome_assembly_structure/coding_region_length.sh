#remove long singleton 'repeat' identified by ROUSFinder, this will generate *.filtered_rep.fas
echo $1
nohup makeblastdb -in temp/$1.fas -dbtype nucl -out mt >/dev/null 2>&1
blastn -db mt -query ~/Downloads/sequence.txt -outfmt 6 -evalue 1e-70 >cds.blast
awk -v OFS='\t' '{if ($9 <= $10) print $2, $9, $10; else print $2, $10, $9}' cds.blast | sort -k1,1 -k2,2n >tmp.bed
bedtools merge -i tmp.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'


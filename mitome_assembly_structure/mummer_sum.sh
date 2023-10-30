echo $1
sed -n '6,$p' $1 | awk -v OFS="\t" '{print "Rgl_1",$1,$2}' | sort -k2,2n >tem.bed
bedtools merge -i tem.bed | awk '{sum += $3 - $2} END {print sum}'

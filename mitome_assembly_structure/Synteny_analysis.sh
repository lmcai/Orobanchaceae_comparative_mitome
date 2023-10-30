#####################################
#MCScan abandoned, this relies on gene synteny to establish genomewise synteny
#makeblastdb -in Aeginetia_indica.fas -dbtype nucl -out temp
#blastn  -task dc-megablast -db temp -query Harveya_capensis.fas -evalue 1e-5 -outfmt 6 -num_threads 4 > AH.blast
#makeblastdb -in Harveya_capensis.fas -dbtype nucl -out temp
#blastn  -task dc-megablast -db temp -query Aeginetia_indica.fas -evalue 1e-5 -outfmt 6 -num_threads 4 >> AH.blast


#####################################
#minimap2-dotplotly-NGenomeSyn
#cross-species genome alignment
./minimap2-2.26_x64-linux/minimap2 -cx asm20 --cs Harveya_capensis.fas Aeginetia_indica.fas >AH.paf

#on MAC, visualize in DotPlotly
#set path to pandoc
#in R use 
#rmarkdown::find_pandoc(cache = FALSE)
#to locate the path of pandoc
export PATH=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/:$PATH
./dotPlotly/pafCoordsDotPlotly.R -i AH.paf -o AH.plot -t -m 50 -q 500 -l -p 5
#-m min aln len
#-q min query len


#####################################
#MUMMER-dotplotly-NGenomeSyn
nucmer -p AH --maxmatch -g 1000 -b 1000 Harveya_capensis.fas Aeginetia_indica.fas 
show-coords -c AH.delta > AH.coords

./dotPlotly/mummerCoordsDotPlotly.R -i AH.coords -o AH -t -m 50 -q 500 -l -p 5

#mummer produces similar result, but include more short scattered aligned blocks

#summarize total length and align size distribution
sed -n '6,$p' RPch.coords| awk -v OFS="\t" '{print "Rgl_1",$1,$2}' >tem.bed
bedtools merge -i tem.bed | awk '{sum += $3 - $2} END {print sum}'

######################################
#plot aln size histogram
#usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

def calculate_n50(numbers):
    sorted_numbers = sorted(numbers, reverse=True)
    total_sum = sum(sorted_numbers)
    fifty_percent_point = total_sum / 2
    running_total = 0
    n50 = None
    for num in sorted_numbers:
        running_total += num
        if running_total >= fifty_percent_point:
            n50 = num
            break
    return n50
    
files=os.listdir('./')
files=[i for i in files if i.endswith('coords')]
for f in files:
	x=open(f).readlines()
	x=[int(l.split()[6]) for l in x[5:]]
	plt.hist(x, bins=30)
	plt.xlabel('Alignment length (bp)')
	plt.ylabel('Frequency')
	plt.title('Aln N50='+str(calculate_n50(x))+'\nAln median='+str(np.median(x)))
	plt.savefig(f.split('.')[0]+'.mt_aln_len.pdf')
	plt.close()
	print(f.split('.')[0], str(calculate_n50(x)), str(np.median(x)))

hemi=os.listdir('hemi')
holo=os.listdir('holo')
for f in hemi:
	x=open('hemi/'+f).readlines()
	x=[int(l.split()[6]) for l in x[5:]]
	x.sort(reverse=True)
	cumulative_sum = [sum(x[:i+1]) for i in range(len(x))]
	plt.plot(list(range(1,len(x)+1)), cumulative_sum, linestyle='-', marker='', color='b',alpha=0.2)



for f in holo:
	x=open('holo/'+f).readlines()
	x=[int(l.split()[6]) for l in x[5:]]
	x.sort(reverse=True)
	cumulative_sum = [sum(x[:i+1]) for i in range(len(x))]
	plt.plot(list(range(1,len(x)+1)), cumulative_sum, linestyle='-', marker='', color='r',alpha=0.2)



#plt.show()
plt.xlabel('Alignment Index')
plt.ylabel('Accumulative length (bp)')
plt.savefig('accumulative_mt_aln_length.pdf')


#############################
#Visualize conserved region in Rehmannia mt

#tem.sh:
```
SP="$1"
sed -n '6,$p' hemi/$1.coords| awk -v OFS="\t" '{print "Rgl_1",$1,$2}' | sed 's/ /\t/g' | sort -k2,2n >tem.bed
bedtools merge -i tem.bed | awk -v col4="$SP" '{print $0,col4}' | sed 's/ /\t/g' >>hemi_synteny.bed
```
sh tem.sh RAgr
sh tem.sh RAse
...

bedtools coverage -d -a mt.bed -b hemi_synteny.bed | awk '{print $4,$5}' >tem.bed
```
#usr/bin/python
lines = open('tem.bed').readlines()

filtered_lines = [lines[0]]  

for i in range(1, len(lines) - 1):
    line = lines[i]
    prev_line = lines[i - 1]
    next_line = lines[i + 1]
    line_parts = line.split()  # Adjust the separator if needed
    if line_parts[1] != prev_line.split()[1] or line_parts[1] != next_line.split()[1]:
        filtered_lines.append(line)


out=open('hemi.coverage2plot.tsv','a')
out.writelines(filtered_lines)
```


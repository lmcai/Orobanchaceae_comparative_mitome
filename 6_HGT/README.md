# HGT

HGTScanner_mt was used to identify homology, build phylogeny, and classify loci.

1. Prepare a bed file for coding regions and MTPT. These region will be masked due to excessive BLAST hits if queried against the NCBI nr/nt database for mitochondrial sequences. Use `gb2gff.py` to convert genbank formated annotation to gff and then to bed format.

See example from `mask.bed`. This file can include all regions that need to be masked across multiple species in the analysis.

2. Run HGTScanner to establish homologous regions for phylogenetic investigation. 

a. Make sure you have `Viridiplantae_mt.taxonomy` and `Viridiplantae_mt.fasta` in your working directory. This is the NCBI mito database. 

b. [optional] Add more mitochondrial genomes (e.g., `oro_mt.fas`) from close relatives to help distinguish VGT vs HGT. Any query sequences that cluster with these close relatives will be classified as VGT. 

c. Run HGTScanner

```
python HGTscanner_mt.py -q [query assembly] -ref oro_mt.fas -o [output prefix] -f Orobanchaceae -b mask.bed
```

d. Output: This will generate a `[prefix]_HGTscanner_supporting_files` folder containing multiple fasta files of homologous regions

3. Sequence alignment and phylogeny

Use the following command to align sequence and infer phylogeny
```
mafft --localpair --maxiterate 1000 --quiet --adjustdirection [prefix].hgt.1.fas | sed 's/_R_//g' > [prefix].hgt.1.aln.fas
iqtree2 -B 1000 -T 4 -redo -s [prefix].hgt.1.aln.fas
```

4. 
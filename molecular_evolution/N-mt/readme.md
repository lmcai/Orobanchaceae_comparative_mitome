N-mt assembly and selection analysis

1. Generate reference of N-mt from public genomes and transcriptomes with PhyloHerb

Orthologs of N-mt from Arabidopsis is obtained from Ceriotti et al (2022). These orthologs are searched against published genomes of Lindenbergia philippensis, Phelipanche, Orobanche, as well as oneKP transcriptomes of Conopholis americana, Epifagus virginiana, and Orobanche fasciculata. 
```
python phyloherb.py -m ortho -i genomes -o ORO_ref -evalue 1e-20 -nuc
```
2. Reformat the reference seqs from Orobanchaceae and retrive these sequences in each genome sequencing reads
```
python phyloherb.py -m assemb -r1 ../cleaned_reads/SRR7688105_1.fastq -r2 ../cleaned_reads/SRR7688105_2.fastq -ref ORO_ref.fas -prefix Cistanchetubulosa -n 8
```

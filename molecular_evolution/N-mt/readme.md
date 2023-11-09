N-mt assembly and selection analysis

1. Generate reference of N-mt from public genomes and transcriptomes with PhyloHerb

Orthologs of N-mt from Arabidopsis is obtained from Ceriotti et al (2022). These orthologs are searched against published genomes of Lindenbergia philippensis, Phelipanche, Orobanche, as well as oneKP transcriptomes of Conopholis americana, Epifagus virginiana, and Orobanche fasciculata. 

Use reciprocal BLAST to identify orthologs with the script `reciprocol_best_hit.py`

2. Extract coding sequences of these transcripts with TransDecoder `transdecoder.sh`. This is necessary because some of the OneKP data and even annotated genomes contain UTRs, which need to be removed for dN/dS.

3. Because some homologs from Orobanchaceae seem to be paralogs. I am adding four more species Glycine, Mumulus, Carica, and Solanum to assist with ortholog assignment. Use similar reciprocal best hit method to get potential orthologs.

4. Reformat the reference seqs from Orobanchaceae and retrive these sequences in each genome sequencing reads
```
python phyloherb.py -m assemb -r1 ../cleaned_reads/SRR7688105_1.fastq -r2 ../cleaned_reads/SRR7688105_2.fastq -ref ORO_ref.fas -prefix Cistanchetubulosa -n 8
```

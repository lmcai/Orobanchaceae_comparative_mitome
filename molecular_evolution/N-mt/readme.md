N-mt assembly and selection analysis

1. Generate reference of N-mt from public genomes and transcriptomes with PhyloHerb

Orthologs of N-mt from Arabidopsis is obtained from Ceriotti et al (2022). These orthologs are searched against published genomes of Lindenbergia philippensis, Phelipanche, Orobanche, as well as oneKP transcriptomes of Conopholis americana, Epifagus virginiana, and Orobanche fasciculata. 

Use reciprocal BLAST to identify orthologs with the script `reciprocol_best_hit.py`

2. Extract coding sequences of these transcripts with TransDecoder `transdecoder.sh`. This is necessary because some of the OneKP data and even annotated genomes contain UTRs, which need to be removed for dN/dS.

3. Because some homologs from Orobanchaceae seem to be paralogs. I am adding four more species Glycine, Mumulus, Carica, and Solanum to assist with ortholog assignment. Use similar reciprocal best hit method to get potential orthologs.

4. Reformat the reference seqs from Orobanchaceae and retrive these sequences in each genome sequencing reads
```
 python ../../PhyloHerb/phyloherb.py -m ortho -i N-mt/ -ref ORO_ref.fas -o test -evalue 1e-20 -nuc
cat: test/AT1G08480.fas: No such file or directory
cat: test/AT1G14450.fas: No such file or directory
cat: test/AT1G24090.fas: No such file or directory
cat: test/AT1G31010.fas: No such file or directory
cat: test/AT1G47720.fas: No such file or directory
cat: test/AT1G71260.fas: No such file or directory
cat: test/AT1G76200.fas: No such file or directory
cat: test/AT2G27730.fas: No such file or directory
cat: test/AT2G46505.fas: No such file or directory
cat: test/AT2G46540.fas: No such file or directory
cat: test/AT3G13226.fas: No such file or directory
cat: test/AT3G47833.fas: No such file or directory
cat: test/AT4G20010.fas: No such file or directory
cat: test/AT4G37830.fas: No such file or directory
cat: test/AT5G13450.fas: No such file or directory
cat: test/AT5G44785.fas: No such file or directory
cat: test/AT5G47890.fas: No such file or directory
cat: test/AT5G50340.fas: No such file or directory
cat: test/AT5G51080.fas: No such file or directory
```

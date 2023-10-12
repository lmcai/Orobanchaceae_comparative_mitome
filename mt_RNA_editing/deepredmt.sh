deepredmt test.fas > seqs.pred
cut -f3,5 seqs.pred | sed 's/C/0/g' | sed 's/E/1/g' > seqs.parsed_pred

#Execute across genes
deepredmt atp1.aln.ordered_NT.fasta >atp1.aln.ordered_NT.pred
deepredmt atp4.aln.ordered_NT.fasta >atp4.aln.ordered_NT.pred
deepredmt atp6.aln.ordered_NT.fasta >atp6.aln.ordered_NT.pred
deepredmt atp8.aln.ordered_NT.fasta >atp8.aln.ordered_NT.pred
deepredmt atp9.aln.ordered_NT.fasta >atp9.aln.ordered_NT.pred
deepredmt ccmB.aln.ordered_NT.fasta >ccmB.aln.ordered_NT.pred
deepredmt ccmC.aln.ordered_NT.fasta >ccmC.aln.ordered_NT.pred
deepredmt ccmFN.aln.ordered_NT.fasta >ccmFN.aln.ordered_NT.pred
deepredmt ccmFc.aln.ordered_NT.fasta >ccmFc.aln.ordered_NT.pred
deepredmt cob.aln.ordered_NT.fasta >cob.aln.ordered_NT.pred
deepredmt cox1.aln.ordered_NT.fasta >cox1.aln.ordered_NT.pred
deepredmt cox2.aln.ordered_NT.fasta >cox2.aln.ordered_NT.pred
deepredmt cox3.aln.ordered_NT.fasta >cox3.aln.ordered_NT.pred
deepredmt matR.aln.ordered_NT.fasta >matR.aln.ordered_NT.pred
deepredmt mttB.aln.ordered_NT.fasta >mttB.aln.ordered_NT.pred
deepredmt nad1_exon1.ordered_NT.fasta >nad1_exon1.ordered_NT.pred
deepredmt nad1_exon23.ordered_NT.fasta >nad1_exon23.ordered_NT.pred
deepredmt nad1_exon4.ordered_NT.fasta >nad1_exon4.ordered_NT.pred
deepredmt nad1_exon5.ordered_NT.fasta >nad1_exon5.ordered_NT.pred
deepredmt nad2_exon12.ordered_NT.fasta >nad2_exon12.ordered_NT.pred
deepredmt nad2_exon345.ordered_NT.fasta >nad2_exon345.ordered_NT.pred
deepredmt nad3.aln.ordered_NT.fasta >nad3.aln.ordered_NT.pred
deepredmt nad4.aln.ordered_NT.fasta >nad4.aln.ordered_NT.pred
deepredmt nad4L.aln.ordered_NT.fasta >nad4L.aln.ordered_NT.pred
deepredmt nad5_exon12.ordered_NT.fasta >nad5_exon12.ordered_NT.pred
deepredmt nad5_exon45.ordered_NT.fasta >nad5_exon45.ordered_NT.pred
deepredmt nad6.aln.ordered_NT.fasta >nad6.aln.ordered_NT.pred
deepredmt nad7.aln.ordered_NT.fasta >nad7.aln.ordered_NT.pred
deepredmt nad9.aln.ordered_NT.fasta >nad9.aln.ordered_NT.pred
deepredmt rpl10.aln.ordered_NT.fasta >rpl10.aln.ordered_NT.pred
deepredmt rpl16.aln.ordered_NT.fasta >rpl16.aln.ordered_NT.pred
deepredmt rpl5.aln.ordered_NT.fasta >rpl5.aln.ordered_NT.pred
deepredmt rps12.aln.ordered_NT.fasta >rps12.aln.ordered_NT.pred
deepredmt rps13.aln.ordered_NT.fasta >rps13.aln.ordered_NT.pred
deepredmt rps14.aln.ordered_NT.fasta >rps14.aln.ordered_NT.pred
deepredmt rps3.aln.ordered_NT.fasta >rps3.aln.ordered_NT.pred
deepredmt rps4.aln.ordered_NT.fasta >rps4.aln.ordered_NT.pred

#Then summarize
awk '{if ($5 >0.9) print FILENAME,$1,$5}' *.deepredmt.pred >deepredmt.all.pred  

#Get sum stat in python
See Mapping_deepredmt_to_alignment.py and summarize_editing_sites.py
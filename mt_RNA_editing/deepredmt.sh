deepredmt test.fas > seqs.pred
cut -f3,5 seqs.pred | sed 's/C/0/g' | sed 's/E/1/g' > seqs.parsed_pred


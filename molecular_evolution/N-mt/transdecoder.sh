TransDecoder.LongOrfs -t $1.aln.fas -m 50
TransDecoder.Predict --single_best_only -t $1.aln.fas

cat $1.aln.fas.transdecoder.cds $1.aln.fas >>$1.tem.fas
mafft $1.tem.fas | awk '/^>/{if (seq) print seq; seq=""; printf $0"\r"; next}{seq = seq $0}END{print seq}' | sort -t$'\t' -k1,1 >$1.aln2.fas

rm $1.tem.fas
rm -r $1.aln.fas.transdecoder_dir
rm -r $1.aln.fas.transdecoder_dir.__checkpoints
rm -r $1.aln.fas.transdecoder_dir.__checkpoints_longorfs

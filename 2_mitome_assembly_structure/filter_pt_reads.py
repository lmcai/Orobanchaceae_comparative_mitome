import pysam

def filter_sam(input_sam, output_sam):
    with pysam.AlignmentFile(input_sam, "r") as samfile_in, \
         pysam.AlignmentFile(output_sam, "wh", template=samfile_in) as samfile_out:
        
        for read in samfile_in:
            if not read.is_unmapped:
                match_length = sum(cigar[1] for cigar in read.cigartuples if cigar[0] == 0)
                read_length = read.query_length
                similarity = match_length / read_length
                
                if similarity < 1.0 or match_length / read_length < 0.98:
                    samfile_out.write(read)

input_sam = "input.sam"
output_sam = "filtered_output.sam"
filter_sam(input_sam, output_sam)

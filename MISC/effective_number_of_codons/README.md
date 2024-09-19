1. Concatenated sequence global (all genes) GC content analysis with `codeml`

Run `codeml codeml.ctl`. Input include `sp43_allmt.phy` and `sp43_allmt.tre`

Then use `parse_codeml_output4codon_usage.py` to extract codon table from the result, as well as GC123 for all genes

2. Individual gene GC content per species

Use neutrality test to examine the correlation between GC12 and GC3.

Use `codon_analysis_fasta_prep.py` to reformat fasta to one file per species and examine the number of consecutive '-' to identify frameshifts.

Then use `codon_GC_F3x4_per_sp.py` to calculate GC12 and GC3 per gene (24 core mt genes per species) for wach species.

Finally, plot GC12 ~ GC3 for each species and conduct Spearman correlation test using `GC123_cor_test_plot.R`
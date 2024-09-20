# Comparative hypothesis testing: Phylogenetic ANOVA and PGLS analysis

## Divergence time estimation in TreePL

The age of six nodes was fixed based on Mortimer et al 2022. TreePL was run three times to prime, cross validate, and thoroughly search for the optimal time tree. Configuration file can be found in `treePL.in.txt`. The output ultrametric tree is `round3.mt37g_42sp.treePL.tre`.

## Hypothesis testing with phylogenetic ANOVA and PGLS

The input data matrix is available in `genomesize_hemiholo_HGT.csv`. The statistical analysis is completed using `PGLS_phylAnova.R`.
# Comparative hypothesis testing: Phylogenetic ANOVA and PGLS analysis

## Divergence time estimation

**1. Orobanchaceae**

The age of six key nodes was fixed based on Mortimer et al 2022. TreePL was run three times to prime, cross validate, and thoroughly search for the optimal time tree. Configuration file can be found in `treePL.in.txt`. The output ultrametric tree is `round3.mt37g_42sp.treePL.tre`.
<img src="./Orobanchaceae_timetree.png" width="570" height="700">

**2. Heterotrophic angiosperms**

A quick and dirty mitochondrial phylogeny of all 19 species was generated using the PhyloHerb pipeline:

```
python ~/Documents/GitHub/PhyloHerb/phyloherb.py -m ortho -i MTPT_other_order -suffix .mt.fasta -mito -o test
```
mtDNA was aligned with MAFFT and concatenated with [catfasta2phyml](https://github.com/nylander/catfasta2phyml), a ML tree was inferred in IQTREE2 with ModelFinder to determine the best substitution model. 
```
mafft atp1.fas >atp1.aln.fas
...

catfasta2phyml.pl -c *.fas > other_heterotrophic_angiosperm.phy 2> partitions.txt
iqtree2 -s other_heterotrophic_angiosperm.phy -p partitions.txt
```

An ultrametric tree was generated in ape with relaxed , root age fixed at 105.4 Ma accoding to Magallón et al. (2015)
```
library(ape)
t=read.tree('other_heterotrophic_angio.roote.tre')
mycalibration <- makeChronosCalib(t, node="root", age.max=105.4)

#correlated model
mytimetree <- chronos(t, lambda = 1, model = "correlated", calibration = mycalibration, control = chronos.control() )
#log-Lik = -1.73748 ; PHIIC = 123.48 

#discrete mode
mytimetree <- chronos(t, lambda = 1, model = "discrete", calibration = mycalibration, control = chronos.control() )
#log-Lik = -1.800698 ; PHIIC = 81.6

#relaxed model
mytimetree <- chronos(t, lambda = 1, model = "relaxed", calibration = mycalibration, control = chronos.control() )
#log-Lik = -1.886355 ; PHIIC = 126.25 
```
The correlated model has the highest log-like and the time tree was selected for PhylANOVA `other_heterotrophic_angio.chronos.tre`.

## Hypothesis testing with phylogenetic ANOVA and PGLS

The input data matrix is available in `hemiholo_traits.csv`. The statistical analysis is completed using `PGLS_phylAnova.R`.

## Coevolution of nucleotide substitution and genomic traits

Coevolution of nucleotide substitution and genomic traits was evaluated using COEVOL v 1.6. It is a Bayesian MCMC program for doing comparative analyses combining molecular data and quantitative traits. 

**1. Prepare input**
a. codon alignment
b. ultrametric tree
c. trait matrix (quantitative traits only)

**2. Running COEVOL**

Combine dS, dN/dS and GC∗ in the same analysis and run two chains
```
coevol -d all_mt_no_ribosome.conc.phy -t round3.mt37g_43sp.treePL.tre -fixtimes -c traits4coevol.txt -dsom -gc orocdsomgc1
```

**3. Check convergence in R**
```#R
library(coda)
x=read.table('orocdsomgc1.log')
mcmc_object <- mcmc(x$V1)
effectiveSize(mcmc_object)
```

**4. Generate covariance matrix**
```readcoevol -x 200 orocdsomgc1
#x is the burn-in
```

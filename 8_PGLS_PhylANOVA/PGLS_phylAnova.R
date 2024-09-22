#####################################################################
#Load library and data
library(ape)
library(nlme)
library(geiger)
data <- read.csv("hemiholo_traits.csv",row.names = 1)
data$Class_code=as.factor(data$Class_code)
phy_tree=read.tree('round3.mt37g_42sp.treePL.tre')
phy_tree$node.label<-NULL
name.check(data,phy_tree)

##########################
#Figure S1: PGLS pt_genome_size ~ gc%
library(caper)
comp.data<-comparative.data(phy_tree, data, names.col="Species_dup", vcv.dim=2)
model_pt<-pgls(pt_size~pt_gc, data=comp.data)

Call:
pgls(formula = pt_size ~ pt_gc, data = comp.data)

Residuals:
     Min       1Q   Median       3Q      Max 
-2.08103 -0.16315  0.01618  0.11618  1.69020 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept) 3.3987e+01 2.3364e+00 14.5471  < 2e-16 ***
pt_gc       2.5111e-05 1.4239e-05  1.7635  0.08565 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.5249 on 39 degrees of freedom
Multiple R-squared: 0.07385,	Adjusted R-squared: 0.0501 
F-statistic:  3.11 on 1 and 39 DF,  p-value: 0.08565 

#######################################################################
#PhyloANOVA of mt_size in hemi vs holo
library(phytools)
data$Class_code=as.factor(data$Class_code)
phylANOVA(phy_tree, data$Class_code, data$mt_size, nsim=1000, posthoc=TRUE, p.adj="holm")

ANOVA table: Phylogenetic ANOVA

Response: y
               Sum Sq      Mean Sq  F value Pr(>F)
x        1.417284e+12 1.417284e+12 3.992117  0.487
Residual 1.384580e+13 3.550206e+11                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -1.998028
1 1.998028  0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.487
1 0.487 1.000
---------

#######################################################################
#Figure 2A: PhyloANOVA of hgt_mtpt in hemi vs holo

data$Class_code=as.factor(data$Class_code)
phylANOVA(phy_tree, data$Class_code, data$mtpt_perc, nsim=1000, posthoc=TRUE, p.adj="holm")

ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq   Mean Sq  F value Pr(>F)
x        377.2589 377.25894 46.12909  0.019
Residual 318.9549   8.17833                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0       1
0  0.00000 6.79184
1 -6.79184 0.00000

Pairwise corrected P-values:
      0     1
0 1.000 0.019
1 0.019 1.000
---------

#######################################################################
#Figure 2A: PhyloANOVA of hgt_mtpt in fully heterotrophic plant in other angiosperm orders

data2 <- read.csv("heterotrophic_angiosperm_mtpt.csv",row.names = 1)
data2$Lifestyle=as.factor(data2$Lifestyle)
phy_tree2=read.tree('other_heterotrophic_angio.chronos.tre')
phy_tree2$node.label<-NULL
name.check(data2,phy_tree2)
phylANOVA(phy_tree2, data2$Lifestyle, data2$mtpt_perc, nsim=1000, posthoc=TRUE, p.adj="holm")

ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq  F value Pr(>F)
x        0.004561 0.004561 7.200591  0.012
Residual 0.010768 0.000633                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
          non-photo     photo
non-photo  0.000000 -2.683392
photo      2.683392  0.000000

Pairwise corrected P-values:
          non-photo photo
non-photo     1.000 0.012
photo         0.012 1.000
---------
##########################
#Figure S8: PGLS hgt ~ mt genome size
#check hgt ~ mt genome size without phylogeny
#see a potential bug solution from https://stackoverflow.com/questions/63138372/error-running-pgls-in-ape-no-covariate-specified

model1<-gls(hgt ~ mt_size, data=data, correlation=corBrownian(1, phy_tree, form = ~Species_dup))
summary(model1)
Generalized least squares fit by REML
  Model: hgt ~ mt_size 
  Data: data 
       AIC      BIC    logLik
  1004.676 1009.666 -499.3378

Correlation Structure: corBrownian
 Formula: ~Species_dup 
 Parameter estimate(s):
numeric(0)

Coefficients:
                Value Std.Error   t-value p-value
(Intercept) -8348.016  49163.98 -0.169799   0.866
mt_size         0.102      0.02  5.313144   0.000

 Correlation: 
        (Intr)
mt_size -0.276

Standardized residuals:
        Min          Q1         Med          Q3         Max 
-0.70967396 -0.13845168 -0.04400526  0.20138983  2.27085583 

Residual standard error: 110939.1 
Degrees of freedom: 41 total; 39 residual

#GLS with caper
library(caper)
data <- read.csv("HGT_genomesize_hemiholo.csv",row.names = 1)
comp.data<-comparative.data(phy_tree, data, names.col="Species_dup", vcv.dim=2)
model2<-pgls(hgt~mt_size, data=comp.data)

Call:
pgls(formula = hgt ~ mt_size, data = comp.data1)

Residuals:
   Min     1Q Median     3Q    Max 
-55300  -5829  -1567   3870  83662 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
               Estimate  Std. Error t value  Pr(>|t|)    
(Intercept) -8.3480e+03  4.9164e+04 -0.1698     0.866    
mt_size      1.0182e-01  1.9164e-02  5.3131 4.643e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 20180 on 39 degrees of freedom
Multiple R-squared: 0.4199,	Adjusted R-squared: 0.405 
F-statistic: 28.23 on 1 and 39 DF,  p-value: 4.643e-06 


###################################################
#Data visualization
#boxplot of consolidated hgt+mtpt content
x=read.csv('hemiholo_traits.csv')
pdf('FigS_hgt_mtpt_content.pdf',width=2,height=4)
boxplot(hgt_mtpt_content ~ Class_code, data=x)
points(jitter(as.numeric(factor(x$Class_code))), x$hgt_mtpt_content, col = "blue", pch = 16)
dev.off()

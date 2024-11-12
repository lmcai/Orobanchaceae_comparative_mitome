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
#PhyloANOVA of syntenic regions compared to Rehmania in hemi vs holo
library(phytools)
data$Class_code=as.factor(data$Class_code)

phylANOVA(phy_tree, data$Class_code, data$synt_w_Reh_size/data$mt_size, nsim=1000, posthoc=TRUE, p.adj="hochberg")
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq  Mean Sq   F value Pr(>F)
x        0.253361 0.253361 21.848527   0.05
Residual 0.463850 0.011596                 

P-value based on simulation.
---------

Pairwise posthoc test using method = "hochberg"

Pairwise t-values:
          0        1
0  0.000000 4.674241
1 -4.674241 0.000000

Pairwise corrected P-values:
     0    1
0 1.00 0.05
1 0.05 1.00
---------
-
#######################################################################
#PhyloANOVA of RNA editing in hemi vs holo
data$Class_code=as.factor(data$Class_code)
phylANOVA(phy_tree, data$Class_code, data$rna_editing, nsim=1000, posthoc=TRUE, p.adj="holm")

ANOVA table: Phylogenetic ANOVA

Response: y
            Sum Sq   Mean Sq  F value Pr(>F)
x         575.6728 575.67283 8.299574  0.289
Residual 2705.1077  69.36173                

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
          0        1
0  0.000000 2.880898
1 -2.880898 0.000000

Pairwise corrected P-values:
      0     1
0 1.000 0.289
1 0.289 1.000
---------

#######################################################################
#Figure 2A: PhyloANOVA of mtpt in hemi vs holo

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
########################################
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

#########################
#HGT percantage between hemi and holo
phylANOVA(phy_tree, data$Class_code, data$hgt_perc, nsim=1000, posthoc=TRUE, p.adj="holm")
ANOVA table: Phylogenetic ANOVA

Response: y
           Sum Sq   Mean Sq   F value Pr(>F)
x        578.3329 578.33286 26.540366   0.048
Residual 871.6276  21.79069                 

P-value based on simulation.
---------

Pairwise posthoc test using method = "holm"

Pairwise t-values:
         0         1
0 0.000000 -5.151734
1 5.151734  0.000000

Pairwise corrected P-values:
     0    1
0 1.00 0.048
1 0.05 1.00
---------

#######################
#Fig.4 holoparasitic species
###########################
#Fig. 4A pt_omega vs pt_omega
library(caper)
phy_tree=read.tree('holo_only.tre')
phy_tree$node.label<-NULL
data <- read.csv("hemiholo_traits.csv",row.names = 1)
data$Class_code=as.factor(data$Class_code)
data=data[data$Species_dup %in% phy_tree$tip.label,]
comp.data<-comparative.data(phy_tree, data, names.col="Species_dup", vcv.dim=2)
model<-pgls(dS_rehmannia~dS_pt_rehmannia, data=comp.data)

Call:
pgls(formula = dS_rehmannia ~ dS_pt_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-8.836e-04 -1.478e-04  6.965e-05  3.033e-04  5.532e-04 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     0.0060486  0.0050566  1.1962   0.3175
dS_pt_rehmannia 0.0368587  0.0203889  1.8078   0.1684

Residual standard error: 0.0006339 on 3 degrees of freedom
Multiple R-squared: 0.5214,	Adjusted R-squared: 0.3618 
F-statistic: 3.268 on 1 and 3 DF,  p-value: 0.1684 

###########################
#Fig. 4B pt_omega vs pt_omega
model<-pgls(dN.dS_pt_rehmannia~dN.dS_rehmannia, data=comp.data)
Call:
pgls(formula = dN.dS_pt_rehmannia ~ dN.dS_rehmannia, data = comp.data)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.006288 -0.003055  0.003540  0.003558  0.003824 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.27588    0.12500  2.2070   0.1144
dN.dS_rehmannia -0.08228    0.29103 -0.2827   0.7958

Residual standard error: 0.005437 on 3 degrees of freedom
Multiple R-squared: 0.02595,	Adjusted R-squared: -0.2987 
F-statistic: 0.07993 on 1 and 3 DF,  p-value: 0.7958 
###############################
#Fig. 4C mt_omega vs mt_dS
model<-pgls(dS_rehmannia~dN.dS_rehmannia, data=comp.data)
Call:
pgls(formula = dS_rehmannia ~ dN.dS_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0005822 -0.0002895  0.0001970  0.0003563  0.0004421 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)  
(Intercept)      0.044032   0.011753  3.7465  0.03320 *
dN.dS_rehmannia -0.070503   0.027364 -2.5765  0.08203 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0005112 on 3 degrees of freedom
Multiple R-squared: 0.6887,	Adjusted R-squared: 0.585 
F-statistic: 6.638 on 1 and 3 DF,  p-value: 0.08203 

###############################
#Fig. 4D pt_omega vs pt_dS
model<-pgls(dS_pt_rehmannia~dN.dS_pt_rehmannia, data=comp.data)
Call:
pgls(formula = dS_pt_rehmannia ~ dN.dS_pt_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0027842 -0.0005895  0.0001376  0.0016717  0.0026085 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)   
(Intercept)        -0.558364   0.062036 -9.0007 0.002895 **
dN.dS_pt_rehmannia  3.228384   0.254687 12.6759 0.001059 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.00243 on 3 degrees of freedom
Multiple R-squared: 0.9817,	Adjusted R-squared: 0.9756 
F-statistic: 160.7 on 1 and 3 DF,  p-value: 0.001059 

#######################
#Fig.4 hemiparasitic species
###############################
phy_tree=read.tree('hemi_only.tre')
phy_tree$node.label<-NULL
data <- read.csv("hemiholo_traits.csv",row.names = 1)
data$Class_code=as.factor(data$Class_code)
data=data[data$Species_dup %in% phy_tree$tip.label,]
comp.data<-comparative.data(phy_tree, data, names.col="Species_dup", vcv.dim=2)

###############################
#Fig. 4A
model<-pgls(dS_rehmannia~dS_pt_rehmannia, data=comp.data)
Call:
pgls(formula = dS_rehmannia ~ dS_pt_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0019256 -0.0011998 -0.0002304  0.0004361  0.0046239 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)     0.006706   0.015297  0.4384   0.6727
dS_pt_rehmannia 0.041263   0.095212  0.4334   0.6762

Residual standard error: 0.002006 on 8 degrees of freedom
Multiple R-squared: 0.02294,	Adjusted R-squared: -0.09919 
F-statistic: 0.1878 on 1 and 8 DF,  p-value: 0.6762 
###############################
#fig. 4B
model<-pgls(dN.dS_pt_rehmannia~dN.dS_rehmannia, data=comp.data)

Call:
pgls(formula = dN.dS_pt_rehmannia ~ dN.dS_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0082348 -0.0020188 -0.0006753  0.0016343  0.0060139 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                 Estimate Std. Error t value  Pr(>|t|)    
(Intercept)      0.224407   0.033170  6.7655 0.0001428 ***
dN.dS_rehmannia -0.108135   0.076472 -1.4140 0.1950631    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.004437 on 8 degrees of freedom
Multiple R-squared:   0.2,	Adjusted R-squared: 0.09996 
F-statistic:     2 on 1 and 8 DF,  p-value: 0.1951 
###############################
#Fig. 4C

model<-pgls(dS_rehmannia~dN.dS_rehmannia, data=comp.data)
Call:
pgls(formula = dS_rehmannia ~ dN.dS_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0022138 -0.0010704 -0.0001652  0.0006482  0.0043851 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)      0.020681   0.014871  1.3907   0.2018
dN.dS_rehmannia -0.019616   0.034285 -0.5722   0.5829

Residual standard error: 0.001989 on 8 degrees of freedom
Multiple R-squared: 0.03931,	Adjusted R-squared: -0.08077 
F-statistic: 0.3274 on 1 and 8 DF,  p-value: 0.5829 

###############################
#Fig. 4D
model<-pgls(dS_pt_rehmannia~dN.dS_pt_rehmannia, data=comp.data)
Call:
pgls(formula = dS_pt_rehmannia ~ dN.dS_pt_rehmannia, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-0.0087878 -0.0034955 -0.0005754  0.0012063  0.0146093 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)  
(Intercept)         0.272105   0.087843  3.0976  0.01472 *
dN.dS_pt_rehmannia -0.681672   0.473087 -1.4409  0.18758  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.006638 on 8 degrees of freedom
Multiple R-squared: 0.206,	Adjusted R-squared: 0.1068 
F-statistic: 2.076 on 1 and 8 DF,  p-value: 0.1876 

###############################
#RNA editing vs ds

phy_tree=read.tree('round3.mt37g_44sp.treePL.tre')
phy_tree$node.label<-NULL
data <- read.csv("hemiholo_traits.csv",row.names = 1)
data$Class_code=as.factor(data$Class_code)
data=data[data$Species_dup %in% phy_tree$tip.label,]
comp.data<-comparative.data(phy_tree, data, names.col="Species_dup", vcv.dim=2)

model<-pgls(dS_rehmannia~rna_editing, data=comp.data)
Call:
pgls(formula = dS_rehmannia ~ rna_editing, data = comp.data)

Residuals:
       Min         1Q     Median         3Q        Max 
-1.711e-03 -8.170e-04  1.609e-05  8.159e-04  1.847e-03 

Branch length transformations:

kappa  [Fix]  : 1.000
lambda [Fix]  : 1.000
delta  [Fix]  : 1.000

Coefficients:
               Estimate  Std. Error t value  Pr(>|t|)    
(Intercept)  0.23153084  0.05074969  4.5622 0.0005331 ***
rna_editing -0.00054927  0.00012767 -4.3022 0.0008597 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.001077 on 13 degrees of freedom
Multiple R-squared: 0.5874,	Adjusted R-squared: 0.5557 
F-statistic: 18.51 on 1 and 13 DF,  p-value: 0.0008597 

###################################################
#Data visualization
#boxplot of consolidated hgt+mtpt content
x=read.csv('hemiholo_traits.csv')
pdf('FigS_hgt_mtpt_content.pdf',width=2,height=4)
boxplot(hgt_mtpt_content ~ Class_code, data=x)
points(jitter(as.numeric(factor(x$Class_code))), x$hgt_mtpt_content, col = "blue", pch = 16)
dev.off()

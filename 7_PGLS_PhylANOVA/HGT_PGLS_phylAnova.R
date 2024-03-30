#boxplot of consolidated hgt+mtpt content
x=read.csv('HGT_genomesize_hemiholo.csv')
pdf('FigS_hgt_mtpt_content.pdf',width=2,height=4)
boxplot(hgt_mtpt_content ~ Class_code, data=x)
points(jitter(as.numeric(factor(x$Class_code))), x$hgt_mtpt_content, col = "blue", pch = 16)
dev.off()

#####################################################################
#PGLS hgt ~ genome size
library(ape)
library(nlme)
library(geiger)
data <- read.csv("HGT_genomesize_hemiholo.csv",row.names = 1)
data$Class_code=as.factor(data$Class_code)
phy_tree=read.tree('round3.mt37g_42sp.treePL.tre')
name.check(data,phy_tree)

#check hgt ~ genome size without phylogeny
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

#####PGLS with caper
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





#######################################################################
#PhyloANOVA of hgt_mtpt in hemi vs holo
library(phytools)

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
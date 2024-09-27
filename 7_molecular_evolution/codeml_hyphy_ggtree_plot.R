#install and load library
library(ape)
library(ggtree)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
library(treeio)
#turns out to be very versatile

#read in codeml output

ml <- read.codeml('~/Dropbox/Parasitic_plant_mitonuclear_coevolution/Molecular_evolution/codeml_free-rate/rst', '~/Dropbox/Parasitic_plant_mitonuclear_coevolution/Molecular_evolution/codeml_free-rate/mt_no_ribo.freerate.codeml.out')

#plot phylogeny with branches colored by dN, dS, or dN/dS

#dN/dS
pdf(file='dN_vs_dS.pdf',width=5,height=8)
ggtree(ml, aes(color=dN_vs_dS)) + 
    scale_color_continuous(name='dN/dS', limits=c(0, 1.5),
    oob=scales::squish, low='dark blue', high='red') +
    theme_tree2(legend.position=c(.9, .2))+geom_tiplab()
dev.off()

#dN
pdf(file='dN.pdf',width=5,height=8)
ggtree(ml, aes(color=dN)) + 
    scale_color_continuous(name='dN',
    oob=scales::squish, low='dark blue', high='red') +
    theme_tree2(legend.position=c(.9, .2))+geom_tiplab()
dev.off()

#dS
pdf(file='dS.pdf',width=5,height=8)
ggtree(ml, aes(color=dS)) + 
    scale_color_continuous(name='dS',
    oob=scales::squish, low='dark blue', high='red') +
    theme_tree2(legend.position=c(.9, .2))+geom_tiplab()
dev.off()

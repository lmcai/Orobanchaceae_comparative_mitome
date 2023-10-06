library(AnaCoDa)
setwd("Dropbox/Parasitic_plant_mitonuclear_coevolution/Molecular_evolution/codon/")
tsv_files <- list.files(pattern = "\\Gene.tsv$")

for (file in tsv_files) {
  # Read the TSV file
  x <- read.table(file, sep = "\t")

  # Perform a correlation test for the first two columns
  correlation_result <- cor.test(x$V4/x$V5,(x$V2+x$V3)/x$V5,method='spearman')
  pdf(paste(file,".GC.pdf",sep=""),width=4,height=4)
  # Create a scatterplot
  plot<-plot(x$V4*300/x$V5,(x$V2+x$V3)*300/(x$V5*2),pch=20,xlab = 'CG3%',ylab = 'GC12%',xlim=c(25,60),ylim=c(38,52))
  
  # Add a regression line
  a=(x$V4*300)/x$V5
  b=(x$V2+x$V3)*300/(x$V5*2)
  abline(lm(b ~ a), col = "blue")

  # Add p-value to the plot
  text(x = 25, y = 51, labels = paste("p-value =", format(correlation_result$p.value, digits = 4)), pos = 4)
  print(plot)
  dev.off()
}

s <- seq(0.25, 0.6, length.out = 100)
> y <- 2 + s + 29 / (s^2 + (1 - s)^2)
> plot(s, y, type = "l", xlab = "s", ylab = "y", main = "Curve y = 2 + x + 29 / (s^2 + (1 - s)^2)")

genome_file <- system.file("extdata", "more_genes.fasta", package = "AnaCoDa")
> test=initializeGenomeObject('/Users/lcai/Dropbox/Parasitic_plant_mitonuclear_coevolution/Molecular_evolution/codon/Aeginetia_indica.codon.fas')
> nc <- getNc(test)
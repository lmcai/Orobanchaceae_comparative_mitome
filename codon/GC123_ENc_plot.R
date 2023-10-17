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





tsv_files <- list.files(pattern = "\\codon.fas$")

for (file in tsv_files) {
	temp=initializeGenomeObject(file)
	nc <- getNc(temp)
	sp=strsplit(file,split="[.]")[[1]][1]
	x=read.table(paste(sp,'.GCperGene.tsv',sep=''),row.name=1)
	for (gene in names(nc)){
		x[gene,'V6']=nc[gene]
	}
	#write codon GC and ENC to file
	write.table(x,paste(sp,'.GCperGene.tsv',sep=''),col.names = c('GC1','GC2','GC3','Seq len','ENC'),sep = '\t')
	#plot to file
	s <- seq(0.25, 0.6, length.out = 100)
	y <- 2 + s + 29 / (s^2 + (1 - s)^2)
	pdf(paste(sp,'.ENC.pdf',sep=''),,width=4,height=4)
	plot<-plot(s, y, type = "l", ylim=c(42,62),xlab='GC3%',ylab='Expected ENc')
	points(x$V4*3/x$V5,x$V6,pch=20)
	print(plot)
	dev.off()
}	
	


#####Plot the combined dataset, color by species or color by gene
x=read.table('combined_GC_ENc.tsv',header = T)
df=data.frame(x[,c('Sp','Gene','ENc')])
df=cbind(df,x$GC3*3/x$Len)
colnames(df) <- c("Sp", "Gene",'Enc','GC3')

curve_function <- function(s) {
  2 + s + 29 / (s^2 + (1 - s)^2)
}

library(ggplot2)

# Create a scatter plot
pdf('ENc_by_species.pdf',width=8,height=6)
ggplot(df, aes(x = GC3, y = Enc, color = Sp)) +
    geom_point() +
    labs(
         x = "GC3",
         y = "ENc") +
    scale_color_discrete(name = "Species")+stat_function(fun = curve_function, color = "red") 
dev.off()

pdf('ENc_by_gene.pdf',width=12,height=6)
ggplot(df, aes(x = GC3, y = Enc, color = Gene)) +
    geom_point() +
    labs(
         x = "GC3",
         y = "ENc") +
    scale_color_discrete(name = "Gene")+stat_function(fun = curve_function, color = "red") 
dev.off()
library(ape)
t=read.tree('MPR_input.tre')

a=x[,2]
names(a)<-x[,1] 
o <- MPR(a, t, "Rehmannia_glutinosa")
lower_est=o[,1]
upper_est=o[,2]
for (i in 3:113){
	a=x[,i]
	names(a)<-x[,1] 
	o <- MPR(a, t, "Rehmannia_glutinosa")
	lower_est=cbind(lower_est,o[,1])
	upper_est=cbind(upper_est,o[,2])
}

write.csv(lower_est,'MPR_lower.csv',col.names = F)
write.csv(upper_est,'MPR_upper.csv',col.names = F)
# plot the curvefile overall for each runs of 20k_40

setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys/Genotype Analysis/090507/LAMARC")
library(Hmisc)

dirs = list.files(pattern="20k_40")
i = sort(dirs)
length(dirs)

windows(width=10, height=8)
par(mfrow=c(3,6), mar=c(4,4,3,0.2), mgp=c(2,.8,0))
j = seq(1, 20000, by=1)
xlims=c(0,2.2)
for(i in 1:length(dirs)){
	k = grep("Pop15", dirs[i])
	head = paste(strsplit(dirs[i], split="")[[1]][1:5], collapse="")
	if(length(k)==0){ # didn't output curvefile for pop15
		curve = read.table(paste(dirs[i],"/curvefile_overall_Theta1.txt", sep=""), skip=3, header=T)
		plot(curve$Ln.Theta1.[j], curve$Likelihood[j], type="l", 
		xlim=xlims, xlab="ln(Theta)", ylab="Likelihood", main=head)
		#text(3, max(curve$Likelihood)*0.7, labels=head)
	} else {
		plot(0,0, col="white", xlim=xlims, xlab="ln(Theta)", ylab="Likelihood", main=head)
	}
}
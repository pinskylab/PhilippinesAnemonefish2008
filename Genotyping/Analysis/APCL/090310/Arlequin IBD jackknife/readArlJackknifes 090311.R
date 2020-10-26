setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008/Analysis/090310/Arlequin IBD jackknife")
source("../../readArlFstBatch 090311.R")
library(vegan)

file = "FSTdist.sum"

cebugeo= matrix(data=c(0,25,50,75,100,150,175,200,225,250,NA,0,25,50,75,125,150,175,200,225, 
NA,NA,0,25,50,100,125,150,175,200,NA,NA,NA,0,25,75,100,125,150,175, 
NA,NA,NA,NA,0,50,75,100,125,150,NA,NA,NA,NA,NA,0,25,50,75,100,NA,NA,NA,NA,NA,NA,0,25,50,75,
NA,NA,NA,NA,NA,NA,NA,0,25,50,NA,NA,NA,NA,NA,NA,NA,NA,0,25,NA,NA,NA,NA,NA,NA,NA,NA,NA,0), nrow=10)

leytegeo = matrix(data=c(0,25,50,100,125,150,175,225,NA,0,25,75,100,125,150,200,
NA,NA,0,50,75,100,125,175,NA,NA,NA,0,25,50,75,125,NA,NA,NA,NA,0,25,50,100,
NA,NA,NA,NA,NA,0,25,75,NA,NA,NA,NA,NA,NA,0,50,NA,NA,NA,NA,NA,NA,NA,0), nrow=8)


fsts = readArlFstBatch(file)

mantelpcebu = numeric()
mantelpleyte = numeric()

mantelrcebu = numeric()
mantelrleyte = numeric()

for(i in 1:13){
	cebu = as.matrix(fsts[[i]][1:10,1:10])
	cebu = cebu/(1-cebu) # linearize fst
	cebu[cebu<0] = 0 # set negatives to zero

	leyte = as.matrix(fsts[[i]][11:18,11:18])
	leyte = leyte/(1-leyte)
	leyte[leyte<0] = 0

	c = mantel(cebugeo, cebu, method="pearson")
	l = mantel(leytegeo, leyte, method="pearson")

	mantelpcebu[i] = c$signif
	mantelpleyte[i] = l$signif

	mantelrcebu[i] = c$statistic
	mantelrleyte[i] = l$statistic
}

windows(width=10, height=8)
#nf=layout(matrix(c(1,15,2,16,3,17,4,18,5,19,6,20,7,21,8,22,9,23,10,24,11,25,12,26,13,0,14,0),4,7))
#layout.show(nf)
par(mfrow=c(3,5))
par(mai=c(0.2,0.3,0.2,0.2))

coefcebu = numeric()
coefleyte = numeric()

for(i in 1:13){
	cebu = as.matrix(fsts[[i]][1:10,1:10])
	cebu = cebu/(1-cebu) # linearize fst
	cebu[cebu<0] = 0 # set negatives to zero

	leyte = as.matrix(fsts[[i]][11:18,11:18])
	leyte = leyte/(1-leyte)
	leyte[leyte<0] = 0

	title = strsplit(strsplit(names(fsts)[i], "no")[[1]][2],".", fixed=T)[[1]][1]
	plot(cebugeo, cebu, main=paste("W/out", title), xaxt="n", yaxt="s", pch=1)
	points(leytegeo, leyte, pch=16)	

	l = lm(cebu[1:100][!is.na(cebu[1:100])] ~ cebugeo[1:100][!is.na(cebugeo[1:100])])
	lines(cebugeo[1:100][!is.na(cebu[1:100])], l$fitted.values, col="blue")
	coefcebu[i] = l$coefficients[2]

	l = lm(leyte[1:100][!is.na(leyte[1:100])] ~ leytegeo[1:100][!is.na(leytegeo[1:100])])
	lines(leytegeo[1:100][!is.na(leyte[1:100])], l$fitted.values, col="orange")
	coefleyte[i] = l$coefficients[2]
}

#par(mai=c(0.3,0.3,0.2,0.2))
#hist(mantelrcebu^2)
#hist(mantelrleyte^2)




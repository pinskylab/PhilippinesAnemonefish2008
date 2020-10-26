# fixed error in 090315 version: no longer set negative Fsts to zero

#setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008/Analysis/090313/Arlequin singleloci")
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/Arlequin singleloci")

source("../../readArlFstBatch 090415.R")
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

for(i in 1:14){
	cebu = as.matrix(fsts[[i]][1:10,1:10])
	cebu = cebu/(1-cebu) # linearize fst

	leyte = as.matrix(fsts[[i]][11:18,11:18])
	leyte = leyte/(1-leyte)

	c = mantel(cebugeo, cebu, method="pearson")
	l = mantel(leytegeo, leyte, method="pearson")

	mantelpcebu[i] = c$signif
	mantelpleyte[i] = l$signif

	mantelrcebu[i] = c$statistic
	mantelrleyte[i] = l$statistic
}

#windows(width=10, height=8)
quartz(width=10, height=8)
#nf=layout(matrix(c(1,15,2,16,3,17,4,18,5,19,6,20,7,21,8,22,9,23,10,24,11,25,12,26,13,0,14,0),4,7))
#layout.show(nf)
par(mfrow=c(3,5))
par(mai=c(0.2,0.3,0.2,0.2))

coefcebu = numeric()
coefleyte = numeric()
r2cebu = numeric()
r2leyte = numeric()

for(i in 1:14){
	cebu = as.matrix(fsts[[i]][1:10,1:10])
	cebu[upper.tri(cebu, diag=T)] = NA
	cebu = cebu/(1-cebu) # linearize fst

	leyte = as.matrix(fsts[[i]][11:18,11:18])
	leyte[upper.tri(leyte, diag=T)] = NA
	leyte = leyte/(1-leyte)

	title = strsplit(strsplit(names(fsts)[i], "CebuLeyte ")[[1]][2],".", fixed=T)[[1]][1]
	plot(cebugeo, cebu, main=title, xaxt="n", yaxt="s", pch=1, ylim = c(-0.03,0.17))
	points(leytegeo, leyte, pch=16)	

	l = lm(cebu[1:100][!is.na(cebu[1:100])] ~ cebugeo[1:100][!is.na(cebugeo[1:100]) & cebugeo[1:100]>0])
	lines(cebugeo[1:100][!is.na(cebu[1:100])], l$fitted.values, col="blue")
	coefcebu[i] = l$coefficients[2]
	r2cebu[i] = summary(l)$r.squared

	l = lm(leyte[1:100][!is.na(leyte[1:100])] ~ leytegeo[1:100][!is.na(leytegeo[1:100]) & leytegeo[1:100]>0])
	lines(leytegeo[1:100][!is.na(leyte[1:100])], l$fitted.values, col="orange")
	coefleyte[i] = l$coefficients[2]
	r2leyte[i] = summary(l)$r.squared
}

#par(mai=c(0.3,0.3,0.2,0.2))
#hist(mantelrcebu^2)
#hist(mantelrleyte^2)

out = data.frame(Leyte_b = coefleyte, Leyte_r2 = r2leyte, Leyte_p = mantelpleyte, Cebu_b = coefcebu, 
Cebu_r2 =r2cebu, Cebu_p = mantelpcebu)
row.names(out) = unlist(strsplit(unlist(strsplit(names(fsts), "CebuLeyte "))[seq(2,28,2)],".", fixed=T))[seq(1,27,2)]


write.csv(out, paste("MantelResults", Sys.Date(),".csv", sep=""))




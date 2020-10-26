# Read in Arlequin jackknifes from 1/12/10
# Calc RMA regression and output to file

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/APCL/090507/Arlequin/")
source("../../../readArlFstBatch 090415.R")
library(vegan)
library(smatr)

file = "jackknife/FSTdist.sum"

# read in geographic distance
geo = read.csv("../Aclarkii_2009-05-14 geo.csv", row.names=1)
cebugeo = as.matrix(geo[1:10,1:10])
leytegeo = as.matrix(geo[11:18, 11:18])

# read in fsts
fsts = readArlFstBatch(file)

num = length(fsts)
mant = data.frame(b = numeric(num), r = numeric(num), p = numeric(num), region=character(num))
mant$region = as.character(mant$region)

for(i in 1:num){
	thesefsts = as.matrix(fsts[[i]])
	thesefsts = thesefsts/(1-thesefsts) # linearize fst

	if(length(grep("Cebu", names(fsts)[i]))>0){ # if Cebu
		m = mantel(cebugeo, thesefsts, method="pearson")
		l = line.cis(y=thesefsts[lower.tri(thesefsts)], x= cebugeo[lower.tri(cebugeo)])
		mant$region[i] = "Cebu"
	}
	if(length(grep("Leyte", names(fsts)[i]))>0){ # if Cebu
		m = mantel(leytegeo, thesefsts, method="pearson")
		l = line.cis(y=thesefsts[lower.tri(thesefsts)], x= leytegeo[lower.tri(leytegeo)])
		mant$region[i] = "Leyte"
	}

	mant$p[i] = m$signif
	mant$r[i] = m$statistic
	mant$b[i] = l$coef[2]
}

row.names(mant) = names(fsts)

write.csv(mant, paste("JackKnifeResults",Sys.Date(),".csv", sep=""))


# diagnostic plots
hist(mant$b[mant$region=="Cebu"], breaks=20)



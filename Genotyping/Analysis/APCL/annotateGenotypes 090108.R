setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008")

surveys = read.csv("GPSSurveys2009-01-08.csv")
locations = read.csv("Collections2008-11-21.csv")
genotypes = read.csv("Genotyping/AFG5,6,7,8 Genotypes wide 090107.csv")
genotypes$ID = as.numeric(t(data.frame(strsplit(as.character(genotypes$Sample.Name), split="L"))[2,]))

# Add 2nd allele for homozygotes
#AC1359
	i = which(!is.na(genotypes$AC1359_1) & is.na(genotypes$AC1359_2))
	genotypes$AC1359_2[i] = genotypes$AC1359_1[i]
#APR_Cf29
	i = which(!is.na(genotypes$APR_Cf29_1) & is.na(genotypes$APR_Cf29_2))
	genotypes$APR_Cf29_2[i] = genotypes$APR_Cf29_1[i]
#APY_65
	i = which(!is.na(genotypes$APY_65_1) & is.na(genotypes$APY_65_2))
	genotypes$APY_65_2[i] = genotypes$APY_65_1[i]
#NNG_012
	i = which(!is.na(genotypes$NNG_012_1) & is.na(genotypes$NNG_012_2))
	genotypes$NNG_012_2[i] = genotypes$NNG_012_1[i]

out = merge(subset(locations, select=c(ID,SurveyNum,lat,long)), genotypes, all.y=T)
out = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), out, all.y=T)

i = order(out$SiteNum, out$SurveyNum, out$ID)

out = out[i,]

# write out with ? in place of NA for CONVERT
write.csv(out, "Genotyping/AFG5,6,7,8 Genotypes wide anno 090107.csv", row.names=F, na="?")
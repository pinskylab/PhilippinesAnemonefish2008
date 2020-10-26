setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")

colstoread = c("NULL", "character", "NULL", "NULL", "character", "NULL", "numeric", "numeric",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL", "NULL")


file1 = read.table("AFG1 Genotypes Table 090305.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file2 = read.table("AFG5,6,7,8 Genotypes Table 090304.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file3 = read.table("AFG11,12,13,14 Genotypes Table trim 090304.csv", header = T, sep =",", colClasses = colstoread[1:34], na.string="?")
file4 = read.table("AFG16 Genotypes Table 090309.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file5 = read.table("AFG17,18,20,21 Genotypes Table 090305.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file6 = read.table("AFG19 Genotypes Table 090208.txt", header = T, sep ="\t", colClasses = colstoread[1:34], na.string="?") # this file is lacking a tab at the end of each line
file7 = read.table("AFG22,25,26,27 Genotypes Table 090305.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file8 = read.table("AFG23,28,29,30 Genotypes Table 090220.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file9 = read.table("AFG24,31,32,33 Genotypes Table 090224.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file10 = read.table("AFG34 Genotypes Table 090309.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file11 = read.table("AFG35-39 Genotypes Table 090303.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file12 = read.table("AFG40-44 Genotypes Table 090305.txt", header=T, sep="\t", colClasses=colstoread, na.string="?")

colnames = c("Sample", "Marker", "A1", "A2")
colnames(file1) = paste(colnames, c("","",".AFG1", ".AFG1"), sep="")
colnames(file2) = paste(colnames, c("","",".AFG5", ".AFG5"), sep="")
colnames(file3) = paste(colnames, c("","",".AFG11", ".AFG11"), sep="")
colnames(file4) = paste(colnames, c("","",".AFG16", ".AFG16"), sep="")
colnames(file5) = paste(colnames, c("","",".AFG17", ".AFG17"), sep="")
colnames(file6) = paste(colnames, c("","",".AFG19", ".AFG19"), sep="")
colnames(file7) = paste(colnames, c("","",".AFG22", ".AFG22"), sep="")
colnames(file8) = paste(colnames, c("","",".AFG23", ".AFG23"), sep="")
colnames(file9) = paste(colnames, c("","",".AFG24", ".AFG24"), sep="")
colnames(file10) = paste(colnames, c("","",".AFG34", ".AFG34"), sep="")
colnames(file11) = paste(colnames, c("","",".AFG35", ".AFG35"), sep="")
colnames(file12) = paste(colnames, c("","",".AFG40", ".AFG40"), sep="")
	
# trim to APCL samples and focal loci
file1 = file1[intersect(grep("APCL",file1$Sample), which(!is.na(file1$A1))),] # remove extra loci from AFG1
file2 = file2[grep("APCL",file2$Sample),]
file3 = file3[intersect(grep("APCL",file3$Sample), which(!is.na(file3$A1))),] # remove extra genotypes from AFG11
file4 = file4[intersect(grep("APCL",file4$Sample), which(!is.na(file4$A1) & file4$Marker != "APY_55" & file4$Marker != "APY_44")),] # remove extra loci from AFG16
file5 = file5[grep("APCL",file5$Sample),]
file6 = file6[intersect(grep("APCL",file6$Sample), which(!is.na(file6$A1))),] # remove extra loci from AFG19
file7 = file7[grep("APCL",file7$Sample),]
file8 = file8[grep("APCL",file8$Sample),]
file9 = file9[grep("APCL",file9$Sample),]
file10 = file10[grep("APCL",file10$Sample),]
file11 = file11[intersect(grep("APCL",file11$Sample), which(!is.na(file11$A1))),] # remove extra loci from AFG35
file12 = file12[intersect(grep("APCL",file12$Sample), which(!is.na(file12$A1))),] # remove extra loci from AFG40

#check for duplicated sample/loci combinations
which(duplicated(paste(file1$Sample, file1$Marker)))
which(duplicated(paste(file2$Sample, file2$Marker)))
which(duplicated(paste(file3$Sample, file3$Marker)))
which(duplicated(paste(file4$Sample, file4$Marker)))
which(duplicated(paste(file5$Sample, file5$Marker)))
which(duplicated(paste(file6$Sample, file6$Marker)))
which(duplicated(paste(file7$Sample, file7$Marker)))
which(duplicated(paste(file8$Sample, file8$Marker)))
which(duplicated(paste(file9$Sample, file9$Marker)))
which(duplicated(paste(file10$Sample, file10$Marker)))
which(duplicated(paste(file11$Sample, file11$Marker)))
which(duplicated(paste(file12$Sample, file12$Marker)))

# do duplicates in file4 and file12 have identical genotypes?
# equalities should return TRUE if all genos are identical
sum(duplicated(paste(file4$Sample, file4$Marker))) == sum(duplicated(paste(file4$Sample, file4$Marker, file4$A1, file4$A2)))
sum(duplicated(paste(file12$Sample, file12$Marker))) == sum(duplicated(paste(file12$Sample, file12$Marker, file12$A1, file12$A2)))

# trim duplicates from file4 (AFG16) and file12 AFG40-44)
file4 = file4[!duplicated(paste(file4$Sample, file4$Marker)),]
file12 = file12[!duplicated(paste(file12$Sample, file12$Marker)),]

# update APY_45 allele calls in AFG16, AFG17,18,20,21, and AFG40-44 to bp
apy45update <- function(x){
	alleles = c(10,12,13,14,15,16,17,18,19,25,27,29)
	bp = c(210,214,216,218,220,222,225,227,229,240,244,246)
	index = match(x,alleles)
	return(bp[index])
}
file4$A1.AFG16[file4$Marker=="APY_45"] = apy45update(file4$A1.AFG16[file4$Marker=="APY_45"])
file4$A2.AFG16[file4$Marker=="APY_45"] = apy45update(file4$A2.AFG16[file4$Marker=="APY_45"])
file5$A1.AFG17[file5$Marker=="APY_45"] = apy45update(file5$A1.AFG17[file5$Marker=="APY_45"])
file5$A2.AFG17[file5$Marker=="APY_45"] = apy45update(file5$A2.AFG17[file5$Marker=="APY_45"])
file12$A1.AFG40[file12$Marker=="APY_45"] = apy45update(file12$A1.AFG40[file12$Marker=="APY_45"])
file12$A2.AFG40[file12$Marker=="APY_45"] = apy45update(file12$A2.AFG40[file12$Marker=="APY_45"])
	


dim(file1)
dim(file2)
dim(file3)
dim(file4)
dim(file5)
dim(file6)
dim(file7)
dim(file8)
dim(file9)
dim(file10)
dim(file11)

geno = merge(file1, file2, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file3, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file4, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file5, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file6, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file7, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file8, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file9, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file10, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file11, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)
geno = merge(geno, file12, all.x=T, all.y=T, by=c("Sample", "Marker"))
dim(geno)

# Remove APCL530 and APCL530.2 due to triple-peak problem
i = geno$Sample != "APCL530" & geno$Sample != "APCL530.2"
geno = geno[i,]

# Fill in homozygotes
alleleones = grep("A1.", names(geno), fixed=T)
alleletwos = grep("A2.", names(geno), fixed=T)
for(i in 1:length(alleleones)){
	j <- !is.na(geno[,alleleones[i]]) & is.na(geno[,alleletwos[i]]) # find homozygotes
	geno[j,alleletwos[i]] <- geno[j,alleleones[i]]
}

# Find consensus genotypes and check for errors
geno$error = NA
geno$numgenos = 0
geno$A1consens = NA
geno$A2consens = NA
alleleones = grep("A1.", names(geno), fixed=T)
alleletwos = grep("A2.", names(geno), fixed=T)

for(i in 1:length(geno$Sample)){
	checkone = alleleones[!is.na(geno[i,alleleones])] # get allele cols that aren't na
	checktwo = alleletwos[!is.na(geno[i,alleletwos])] # get allele cols that aren't na

	if(length(checkone) != length(checktwo)){ print(paste("A1 and A2 not equal length", i))}

	if(length(checkone)>0 & length(checktwo)>0){
		matchone = all(geno[i,checkone] == geno[i,checkone[1]]) # check for equality
		matchtwo = all(geno[i,checktwo] == geno[i,checktwo[1]]) # check for equality

		geno$numgenos[i] = length(checkone)

		if(matchone & matchtwo){ 
			if(length(checkone)>1){	
				geno$error[i] = FALSE
			}
			geno$A1consens[i] = geno[i,checkone[1]] 
			geno$A2consens[i] = geno[i,checktwo[1]]
		}
		if(!(matchone & matchtwo)){ 
			geno$error[i] = TRUE
		}
	}
	
	if((i %% 500) == 0){ print(i)}
}

i=table(geno$error, geno$Marker)
i
round((i[2,]/(i[1,]+i[2,]))*100,1)


# Write out errors
i = which(geno$error==T)
write.csv(geno[i,], paste("Genotyping_errors_", Sys.Date(), ".csv", sep=""))

# Calc % complete
i = table(!is.na(geno$A1consens), geno$Marker)
i
round((i[2,]/(i[1,]+i[2,]))*100,1)


# Calc # genos missing per individual
i = table(!is.na(geno$A1consens), geno$Sample)
max(i[1,])
j = matrix(data= c(0,1,2,3,4,5,6), nrow=3, ncol=7, byrow=T)
j[2,] = hist(i[1,], plot=F, breaks = c(-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5,6.5))$counts
j[3,] = round(j[2,]/378*100)
rownames(j) = c("Num Missing", "Num Samples", "% Samples")
print(j)


# Fix genotyping errors where calls were fixable
geno$A1consens[geno$Sample == "APCL166" & geno$Marker == "ACH_B9"] = 279
geno$A2consens[geno$Sample == "APCL166" & geno$Marker == "ACH_B9"] = 313

geno$A1consens[geno$Sample == "APCL182" & geno$Marker == "ACH_B9"] = 313
geno$A2consens[geno$Sample == "APCL182" & geno$Marker == "ACH_B9"] = 365

geno$A1consens[geno$Sample == "APCL125" & geno$Marker == "ACH_D1"] = 334
geno$A2consens[geno$Sample == "APCL125" & geno$Marker == "ACH_D1"] = 334

geno$A1consens[geno$Sample == "APCL274" & geno$Marker == "APR_Cf42"] = 333
geno$A2consens[geno$Sample == "APCL274" & geno$Marker == "APR_Cf42"] = 501

geno$A1consens[geno$Sample == "APCL406" & geno$Marker == "APR_Cf42"] = 389
geno$A2consens[geno$Sample == "APCL406" & geno$Marker == "APR_Cf42"] = 509


# Check loci for regular alleles (no alleles 1bp apart, should match bins in GeneMapper)
loci = unique(geno$Marker)
for(i in loci){
	j=sort(unique(c(geno$A1consens[geno$Marker == i], geno$A2consens[geno$Marker== i])))
	cat(c(i,j))
	cat("\n")
}

# Write out whole table
write.csv(geno, paste("Aclarkii_genotypes_", Sys.Date(), ".csv", sep=""))



# Format for GenAlEx 6
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,lat,long)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

	# sort by region and sample ID
	i = order(genowide$Region, genowide$SiteNum, genowide$Sample)
	genowide = genowide[i,]

	# remove ACH_A7, ACH_B9, APR_Cf42 because of high error and low completeness (3/9/09)
	widenames = names(genowide)
	i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
	genowide = genowide[,-i]
	
	
	# Write to file
	genalexfile = file(paste("Aclarkii_GenAlEx_", Sys.Date(), ".csv", sep=""))
	open(genalexfile, "w")

	# Header: num loci, num samples, num pops, # indiv in each pop, # indivs in each region
	outline = c(13, dim(genowide)[1], length(unique(genowide$SiteNum)), table(genowide$SiteNum), table(genowide$Region))
	cat(outline, file=genalexfile, sep=",", append=F)
	cat("\n", file=genalexfile, append=T)
	
	outline = c("Aclarkii 3 regions", "", "", names(table(genowide$SiteNum)), names(table(genowide$Region)))
	cat(outline, file=genalexfile, sep=",", append=T)
	cat("\n", file=genalexfile, append=T)

	i = grep("A1", names(genowide))
	j = gsub("A1.", "", names(genowide)[i], fixed=T)
	j = paste(j, collapse=",,")
	outline = paste("Sample no.,Pop,",j, ",,,Lat,Long",sep="") 
	cat(outline, file=genalexfile, append=T)
	cat("\n", file=genalexfile, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[,c(grep("Sample",widenames),grep("SiteNum",widenames),grep("A[[:digit:]]",widenames),grep("blank",widenames),grep("lat",widenames), grep("long",widenames))]
	write.table(outfile, file=genalexfile, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(genalexfile)	
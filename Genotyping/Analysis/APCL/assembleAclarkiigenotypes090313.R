setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")



### Assemble the genotypes
colstoread = c("NULL", "character", "NULL", "NULL", "character", "NULL", "numeric", "numeric",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL",  "NULL", "NULL")


file1 = read.table("AFG1 Genotypes Table 090305.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file2 = read.table("AFG5,6,7,8 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file3 = read.table("AFG11,12,13,14 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file4 = read.table("AFG16 Genotypes Table 090309.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file5 = read.table("AFG17,18,20,21 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file6 = read.table("AFG19 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file7 = read.table("AFG22,25,26,27 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file8 = read.table("AFG23,28,29,30 Genotypes Table 090313.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file9 = read.table("AFG24,31,32,33 Genotypes Table 090313.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file10 = read.table("AFG34 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file11 = read.table("AFG35-39 Genotypes Table 090312.txt", header = T, sep ="\t", colClasses = colstoread, na.string="?")
file12 = read.table("AFG40-44 Genotypes Table 090312.txt", header=T, sep="\t", colClasses=colstoread, na.string="?")

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
i = rbind(i,round((i[2,]/(i[1,]+i[2,]))*100,1))
rownames(i) = c("No Errors", "Errors", "Error %")
t(i)

# Write out errors
i = which(geno$error==T)
write.csv(geno[i,], paste("Genotyping_errors_", Sys.Date(), ".csv", sep=""))

# Calc % complete
i = table(!is.na(geno$A1consens), geno$Marker)
i = rbind(i,round((i[2,]/(i[1,]+i[2,]))*100,1))
rownames(i) = c("Missing", "Complete", "% Complete")
t(i)


# Calc # genos missing per individual
i = table(!is.na(geno$A1consens), geno$Sample)
max(i[1,])
j = matrix(data= c(0,1,2,3,4,5,6), nrow=3, ncol=7, byrow=T)
j[2,] = hist(i[1,], plot=F, breaks = c(-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5,6.5))$counts
j[3,] = round(j[2,]/378*100)
rownames(j) = c("Num Missing", "Num Samples", "% Samples")
print(j)


# Calc # genos missing per individual (w/out APR_Cf42 and ACH_B9)
temp = geno
dim(temp)
i = c(grep("ACH_B9", temp$Marker), grep("APR_Cf42", temp$Marker))
temp = temp[-i,]
dim(temp)
i = table(!is.na(temp$A1consens), temp$Sample)
max(i[1,])
j = matrix(data= c(0,1,2,3,4,5,6), nrow=3, ncol=7, byrow=T)
j[2,] = hist(i[1,], plot=F, breaks = c(-0.1, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5,6.5))$counts
j[3,] = round(j[2,]/378*100)
rownames(j) = c("Num Missing", "Num Samples", "% Samples")
print(j)


# Fix genotyping errors where calls were fixable
# (none left to fix)

# Check loci for regular alleles (no alleles 1bp apart, should match bins in GeneMapper)
loci = unique(geno$Marker)
for(i in loci){
	j=sort(unique(c(geno$A1consens[geno$Marker == i], geno$A2consens[geno$Marker== i])))
	cat(c(i,j))
	cat("\n")
}


###################################################
# Write out whole table
write.csv(geno, paste("Aclarkii_genotypes_", Sys.Date(), ".csv", sep=""))
###################################################

########################
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

	# remove ACH_B9 and APR_Cf42 because of high error and low completeness (3/9/09)
	widenames = names(genowide)
	i = c(grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
	genowide = genowide[,-i]
	
	
	# Write to file
	genalexfile = file(paste("Aclarkii_GenAlEx_", Sys.Date(), ".csv", sep=""))
	open(genalexfile, "w")

	# Header: num loci, num samples, num pops, # indiv in each pop, # regions, # indivs in each region
	outline = c(14, dim(genowide)[1], length(unique(genowide$SiteNum)), table(genowide$SiteNum), length(unique(genowide$Region)), table(genowide$Region))
	cat(outline, file=genalexfile, sep=",", append=F)
	cat("\n", file=genalexfile, append=T)
	
	outline = c("Aclarkii 3 regions", "", "", names(table(genowide$SiteNum)), "", names(table(genowide$Region)))
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
	
	
	
###############################################
## Divide by size classes, output to GenAlEx

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)

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
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)


hist(genowide$Size)
quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]

sum(first)
sum(second)
sum(third)
sum(fourth)

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4

# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

	# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
	widenames = names(genowide)
	i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
	genowide = genowide[,-i]
	
	
	# Write to file
	genalexfile = file(paste("Aclarkii_sizes_", Sys.Date(), "_GX.csv", sep=""))
	open(genalexfile, "w")

	# Header: num loci, num samples, num pops, # indiv in each pop, # regions, # indivs in each region
	outline = c(13, dim(genowide)[1], length(unique(genowide$sizeclass)), table(genowide$sizeclass))
	cat(outline, file=genalexfile, sep=",", append=F)
	cat("\n", file=genalexfile, append=T)
	
	outline = c("Aclarkii sizeclasses", "", "", names(table(genowide$sizeclass)))
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
	outfile = genowide[,c(grep("Sample",widenames),grep("sizeclass",widenames),grep("A[[:digit:]]",widenames),grep("blank",widenames),grep("lat",widenames), grep("long",widenames))]
	write.table(outfile, file=genalexfile, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(genalexfile)
	
	
	
###############################################
## Output to MLNE2: All samples (lower quantile vs. larger quantile, all populations combined)

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	
	
	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]

sum(first)
sum(second)
sum(third)
sum(fourth)

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4

# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]
	

# Write to file
	mlnefile = file(paste("Aclarkii_", Sys.Date(), "_MLNE.csv", sep=""))
	open(mlnefile, "w")

	# Header: 0 (closed pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(0, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(35000, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(3, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(0, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# num alleles per locus (only for adults and juvs, sizeclasses 4 or 1)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = genowide$sizeclass == 1 | genowide$sizeclass == 4
	for(i in cols){ # for each locus
		j = union(genowide[a,i], genowide[a,i+1])
		j = j[!is.na(j)]
		alleles[match(i,cols)] = list(sort(j))
		outline = c(outline, length(j))
	}
	cat(paste(outline, collapse=","),file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
			
	# num samples from focal pop
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# generations when samples were taken
	cat("0,1", file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	
	# For adults: Numbers of copies of each allele, at each locus, from each sample of the focal population
	cols = grep("A1.", names(genowide))
	a = genowide$sizeclass == 4
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	# For juvs
	cols = grep("A1.", names(genowide))
	a = genowide$sizeclass == 1
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	close(mlnefile)	


###############################################
## Output to MLNE2: All samples (lower quantile vs. larger quantile, two pops as focal, others as source)

focal = 22
focal2 = 23
focal3 = 24
focal4 = 25

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	
	
	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]

sum(first)
sum(second)
sum(third)
sum(fourth)

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4

# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]
	

# Write to file
	mlnefile = file(paste("Aclarkii_", Sys.Date(), "_MLNEfocal", focal, focal2, focal3, focal4, ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 1 (open pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(1, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(35000, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(3, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(0, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# num alleles per locus (only for adults and juvs, sizeclasses 4 or 1)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = genowide$sizeclass == 1 | genowide$sizeclass == 4
	for(i in cols){ # for each locus
		j = union(genowide[a,i], genowide[a,i+1])
		j = j[!is.na(j)]
		alleles[match(i,cols)] = list(sort(j))
		outline = c(outline, length(j))
	}
	cat(paste(outline, collapse=","),file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
			
	# num samples from focal pop
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# generations when samples were taken
	cat("0,1", file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	
	# For adults in focal pops:
	# Numbers of copies of each allele, at each locus, from each sample of the focal population
	cols = grep("A1.", names(genowide))
	a = genowide$sizeclass == 4 & (genowide$SiteNum == focal | genowide$SiteNum == focal2 | genowide$SiteNum == focal3 | genowide$SiteNum == focal4)
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	# For juvs in focal pop
	cols = grep("A1.", names(genowide))
	a = genowide$sizeclass == 1 & (genowide$SiteNum == focal | genowide$SiteNum == focal2 | genowide$SiteNum == focal3 | genowide$SiteNum == focal4)
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}

	# For source pop (use all samples not from focal)
	cols = grep("A1.", names(genowide))
	a = genowide$SiteNum != focal & genowide$SiteNum != focal2 & genowide$SiteNum != focal3 & genowide$SiteNum != focal4
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	# number of starting points
	cat("1", file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
		
	
	close(mlnefile)	



###############################################
## Output for Geneland

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	

gnlndfile =(paste("Aclarkii_", Sys.Date(), "_Gnlndgen.csv", sep=""))
geofile =(paste("Aclarkii_", Sys.Date(), "_Gnlndgeo.csv", sep=""))

# Write out genetic data
write.table(subset(genowide, select=c(-Sample)), gnlndfile, sep=" ", row.names=F, col.names=F)
	
	# add lat/long info
	i = genowide$Sample
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(genowide, subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), all.x=T, by="Sample", sort=F)
	dim(genowide)
	
	identical(i, genowide$Sample) # must be TRUE: same order of samples
	
	
# Write out geographic data in Lambert coordinates
## load packages 
require(mapproj) 
require(maps) 
## check 
plot(genowide$long, genowide$lat,type="n",xlab="Lon",ylab="Lat",asp=1) 
points(genowide$long, genowide$lat,col=2) 
map(resol=0,add=TRUE) 
## convert (Lon,Lat) coordinates into Lambert 
mapproj.res <- mapproject(x=genowide$long, y=genowide$lat, projection="lambert", 
param=c(min(genowide$lat),max(genowide$lat))) 
## save planar coordinates as a two-column matrix 
coord.lamb <- cbind(mapproj.res$x,mapproj.res$y) 
colnames(coord.lamb) = c("X", "Y")
write.table(coord.lamb, geofile, row.names=F, col.names=F)



###############################################
## Output for structure

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	


	# Write to file
	structfile =file(paste("Aclarkii_", Sys.Date(), "_struct.csv", sep=""))
	open(structfile, "w")

	# Header: loci names
	i = grep("A1", names(genowide))
	j = gsub("A1.", "", names(genowide)[i], fixed=T)
	j = paste(j, collapse="  ")
	j = paste(" ",j,sep="")
	cat(j, file= structfile, sep=" ", append=F)
	cat("\n", file= structfile, append=T)

	# Data: Sample, cols of genotype data
	widenames = names(genowide)
	outfile = genowide[,c(grep("Sample",widenames),grep("A[[:digit:]]",widenames))]
	write.table(outfile, file= structfile, append=T, quote=F,row.names=F, col.names=F, sep=" ", na="-99")
	
	
	close(structfile)

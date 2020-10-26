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
## Output to MLNE2: All samples (lower quantile vs. larger quantile, four pops as focal, others as source)
focalpops = data.frame(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27))


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
fourth_inclusive = genowide$Size >= quants[3]

#sum(first)
#sum(second)
#sum(third)
#sum(fourth)

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4

genowide$sizeclassinclusive = NA
genowide$sizeclassinclusive[first] = 1
genowide$sizeclassinclusive[second] = 2
genowide$sizeclassinclusive[third] = 3
genowide$sizeclassinclusive[fourth_inclusive] = 4

# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

# How many adults and juvs in each population group? Use sizeclassinclusive so that we get enough adults
for(rep in 1:(dim(focalpops)[2])){
	
	focal = focalpops[1,rep]
	focal2 = focalpops[2,rep]
	focal3 = focalpops[3,rep]
	focal4 = focalpops[4,rep]

	a = genowide$sizeclassinclusive == 4 & (genowide$SiteNum == focal | genowide$SiteNum == focal2 | genowide$SiteNum == focal3 | genowide$SiteNum == focal4)

	b = genowide$sizeclass == 1 & (genowide$SiteNum == focal | genowide$SiteNum == focal2 | genowide$SiteNum == focal3 | genowide$SiteNum == focal4)
	print(paste("Rep ", rep, ": Num Adults: ", sum(a), "   Num Juvs: ", sum(b), sep=""))

}

# Write all groups to files
for(rep in 1:(dim(focalpops)[2])){
	
	focal = focalpops[1,rep]
	focal2 = focalpops[2,rep]
	focal3 = focalpops[3,rep]
	focal4 = focalpops[4,rep]
	
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
	a = genowide$sizeclass == 1 | genowide$sizeclassinclusive == 4
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
	a = genowide$sizeclassinclusive == 4 & (genowide$SiteNum == focal | genowide$SiteNum == focal2 | genowide$SiteNum == focal3 | genowide$SiteNum == focal4)
	print(paste("Rep ", rep, ": Num Adults: ", sum(a), sep=""))
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
}




###############################################
## Output to MLNE2: All samples (TopTwo on anemone >= 8cm vs. lower quartile, some pops as focal, others as source)

# sets of 4
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27))
prefix="TopTwo"

# sets of 4, 3, 2, and 1 (for testing only)
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27), m = c(7,8,9), n = c(8,9,10), o=c(9,10,11), p=c(10,11,13), q=c(11,13,14), r=c(13,14,15), s=c(14,15,16), t=c(15,16,17), u=c(18,19,20), v=c(19,20,22), w=c(20,22,23), x=c(22,23,24), y=c(23,24,25), z=c(24,25,27), aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27), qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)

# sets of 1,2, and 3, all using pop 25 (determined from sets above)
focalpops = list(a = 25, b = c(24,25), c=c(23,24,25))
prefix="TopTwo"

# sets of 2
focalpops = list(aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27))
prefix="TopTwo"

# sets of 1
focalpops = list(qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)
prefix="TopTwo"


NeMax = 10000

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
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
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
for(rep in 1:length(focalpops)){	
	focal = focalpops[[rep]]	
	sites = matchall(focal, genowide$SiteNum)
	a = intersect(which(genowide$TopTwo & genowide$Size >= 8), sites)
	b = intersect(which(genowide$sizeclass == 1), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"

# Write all groups to files
for(rep in 1:length(focalpops)){
	
	focal = focalpops[[rep]]

	# Write to file
	mlnefile = file(paste("MLNE/Aclarkii_", Sys.Date(), "_MLNEfocal", prefix, paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 1 (open pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(1, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(0, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# num alleles per locus (only for adults and juvs)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = genowide$AdJuv == "Ad" | genowide$AdJuv == "Juv"
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
	
	sites = matchall(focal, genowide$SiteNum)

	# For adults in focal pops:
	# Numbers of copies of each allele, at each locus, from each sample of the focal population
	cols = grep("A1.", names(genowide))
	a = intersect(which(genowide$AdJuv == "Ad"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), sep=""))
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
	a = intersect(which(genowide$AdJuv == "Juv"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Juvs: ", length(a), sep=""))
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
	a = setdiff(1:length(genowide$SiteNum), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Source: ", length(a), sep=""))
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
}



###############################################
## Output to MLNE2: Focal pops as closed (TopTwo on anemone >= 8cm vs. lower quartile)

# sets of 4
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27))
prefix="TopTwo"

# sets of 4, 3, 2, and 1 (for testing only)
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27), m = c(7,8,9), n = c(8,9,10), o=c(9,10,11), p=c(10,11,13), q=c(11,13,14), r=c(13,14,15), s=c(14,15,16), t=c(15,16,17), u=c(18,19,20), v=c(19,20,22), w=c(20,22,23), x=c(22,23,24), y=c(23,24,25), z=c(24,25,27), aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27), qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)

# sets of 1,2, and 3, all using pop 25 (determined from sets above)
focalpops = list(a = 25, b = c(24,25), c=c(23,24,25))
prefix="TopTwo"

# sets of 2
focalpops = list(aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27))
prefix="TopTwo"

# sets of 1
focalpops = list(qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)
prefix="TopTwo"


NeMax = 10000

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
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
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
for(rep in 1:length(focalpops)){	
	focal = focalpops[[rep]]	
	sites = matchall(focal, genowide$SiteNum)
	a = intersect(which(genowide$TopTwo & genowide$Size >= 8), sites)
	b = intersect(which(genowide$sizeclass == 1), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"

# Write all groups to files
for(rep in 1:length(focalpops)){
	
	focal = focalpops[[rep]]

	# Write to file
	mlnefile = file(paste("MLNE/Aclarkii_", Sys.Date(), "_MLNEclosed", prefix, paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 0 (closed pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(0, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(1, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	sites = matchall(focal, genowide$SiteNum)

	# num alleles per locus (only for adults and juvs in focal pop)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = intersect(which(genowide$AdJuv == "Ad" | genowide$AdJuv == "Juv"), sites)
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
	a = intersect(which(genowide$AdJuv == "Ad"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), sep=""))
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
	a = intersect(which(genowide$AdJuv == "Juv"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Juvs: ", length(a), sep=""))
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
}




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


###############################################
## Output for Alleles in Space (Miller 2005 J Hered)

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
	

aisfile =(paste("Aclarkii_", Sys.Date(), "_AISgen.csv", sep=""))
geofile =(paste("Aclarkii_", Sys.Date(), "_AISgeo.csv", sep=""))

# Write out genetic data
write.table(13,aisfile, sep="", row.names=F, col.names=F, append=F) # num loci

# collapse locus alleles with \
genowide[is.na(genowide)] = 0
cols = grep("A1.", names(genowide))
out = as.character(genowide$Sample)
out = cbind(out, paste(genowide[,cols[1]], genowide[,(cols[1]+1)], sep="\\"))
for(i in cols[2:length(cols)]){
	out = cbind(out, paste(genowide[,i], genowide[,(i+1)], sep="\\"))
}
write.table(out, aisfile, sep=", ", row.names=F, col.names=F, append=T, quote=F)
write.table(";", aisfile, sep=", ", row.names=F, col.names=F, append=T, quote=F)
	
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
offset=1 # add an arbitrary offset so that all coords are positive
coord.lamb <- cbind(as.character(genowide$Sample), (mapproj.res$x+offset),(mapproj.res$y+offset)) 
write.table(coord.lamb, geofile, row.names=F, col.names=F, sep=", ", quote=F, append=F)
write.table(";", geofile, sep=", ", row.names=F, col.names=F, append=T, quote=F)




###############################################
## LDNE: Output juvs and adults separately: TopTwo on anemone >= 8cm vs. lower quartile)
## Uses Genepop format

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
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
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"


# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
groups =c("Ad", "Juv")
pops = sort(unique(genowide$SiteNum))
for(j in 1:length(pops)){	
	a = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Ad")
	b = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Juv")
	print(paste("Pops ", pops[j], ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

# Write all groups to files
cols = grep("A1.", names(genowide))
loci = unlist(strsplit(names(genowide)[cols], "A1.", fixed=T))[seq(2,26, by=2)]
for(i in 1:2){
	
	# Write to file
	ldnename = paste("Aclarkii_", Sys.Date(), "_LDNE", groups[i], ".gen", sep="")
	ldnefile = file(ldnename)
	open(ldnefile, "w")

	# Header
	cat(paste("Title line: Aclarkii all pops, ", groups[i], "s only", sep=""), file= ldnefile, sep=",", append=F)
	cat("\n", file= ldnefile, append=T)
	cat(paste(loci, collapse="\n"), file= ldnefile, append=T)
	cat("\n", file= ldnefile, append=T) 

	for(j in 1:length(pops)){
		k = which(genowide$SiteNum==pops[j] & genowide$AdJuv == groups[i])
		if(length(k)>0){
			cat("Pop", file= ldnefile, sep="", append=T)
			cat("\n", file= ldnefile, append=T)
		
			# collapse locus alleles with ""
			out = paste("     ", rep(pops[j], length(k)), ",", sep="") # name each indiv according to population name
			out = cbind(out, paste(formatC(genowide[k,cols[1]], width=3, flag="0"), formatC(genowide[k,(cols[1]+1)], width=3, flag="0"), sep=""))
			for(c in cols[2:length(cols)]){
				out = cbind(out, paste(formatC(genowide[k,c], width=3, flag="0"), formatC(genowide[k,(c+1)], width=3, flag="0"), sep=""))
			}
			out[out==" NA NA"] = "000000"
			for(c in 1:dim(out)[1]){
				cat(paste(out[c,], collapse=" "), file=ldnefile, append=T)
				cat("\n", file= ldnefile, append=T)
			}	
		}
	}
	close(ldnefile)	
}




###############################################
## Output to ONeSAMP
## Uses a modified Genepop format

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
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
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

# Convert to perfect repeat loci
cols = grep("A1.", names(genowide))
loci = unlist(strsplit(names(genowide)[cols], "A1.", fixed=T))[seq(2,26, by=2)]
for(i in 1:length(loci)){
	j = sort(unique(c(genowide[,cols[i]], genowide[,cols[i]+1])))
	print(paste(loci[i], paste(j, collapse=" ")))
}
	# Rule: shift alleles to 1 lower if needed
	# Perfect di: AC1359, AC1578, AC915, APR_Cf29, APR_Cf8, APY_65, NNG_004, NNG_012
	# Perfect tetra: ACH_D1 
	# Imperfect tetra: ACH_B6: allele 129 -> 127
		i = which(genowide$A1.ACH_B6==129)
		genowide$A1.ACH_B6[i] = 127
		i = which(genowide$A2.ACH_B6==129)
		genowide$A2.ACH_B6[i] = 127		
	# Imperfect di: ACH_C1: 101->100, 105->104, 113->112, 121->120
		i = which(genowide$A1.ACH_C1>100)
		genowide$A1.ACH_C1[i] = genowide$A1.ACH_C1[i]-1
		i = which(genowide$A2.ACH_C1>100)
		genowide$A2.ACH_C1[i] = genowide$A2.ACH_C1[i]-1
	# Imperfect di: APY_45: 225->224, 227->226, 229->228
		i = which(genowide$A1.APY_45>224 & genowide$A1.APY_45<230)
		genowide$A1.APY_45[i] = genowide$A1.APY_45[i]-1
		i = which(genowide$A2.APY_45>224 & genowide$A2.APY_45<230)
		genowide$A2.APY_45[i] = genowide$A2.APY_45[i]-1
	# Imperfect di: NNG_028: 288->287, 290->289, 292->291, 294->293, 296->295, 298->297, 300->299
		i = which(genowide$A1.NNG_028<301)
		genowide$A1.NNG_028[i] = genowide$A1.NNG_028[i]-1
		i = which(genowide$A2.NNG_028<301)
		genowide$A2.NNG_028[i] = genowide$A2.NNG_028[i]-1

# get repeat motif and check for perfect repeats again
repeats = numeric(0)
for(i in 1:length(loci)){
	j = sort(unique(c(genowide[,cols[i]], genowide[,cols[i]+1])))
	temp = numeric(0)
	for(k in 2:length(j)){
		temp = c(temp, j[k] - j[k-1])
	}
	min = min(temp)
	mods = unique(temp %% min)
	if(mods==0){
		repeats = c(repeats,min)
	} else {
		print("locus isn't perfect!")	
	}
}	 
	
# Pop to write out
pop = 22

# Write to file
	name = paste("Aclarkii_", Sys.Date(), "_ONeSAMP", pop, ".gen", sep="")
	file = file(name)
	open(file, "w")

	# Header
	cat(paste("Title line: Aclarkii all pops for ONeSAMP", sep=""), file= file, sep=",", append=F)
	cat("\n", file= file, append=T)
	cat(paste(loci, ", ", repeats, sep="", collapse="\n"), file= file, append=T)
	cat("\n", file= file, append=T) 

		k = which(genowide$SiteNum==pop)
		if(length(k)>0){
			cat("Pop", file= file, sep="", append=T)
			cat("\n", file=file, append=T)
		
			# collapse locus alleles with ""
			out = paste("     ", rep(pop, length(k)), ",", sep="") # name each indiv according to population name
			out = cbind(out, paste(formatC(genowide[k,cols[1]], width=3, flag="0"), formatC(genowide[k,(cols[1]+1)], width=3, flag="0"), sep=""))
			for(c in cols[2:length(cols)]){
				out = cbind(out, paste(formatC(genowide[k,c], width=3, flag="0"), formatC(genowide[k,(c+1)], width=3, flag="0"), sep=""))
			}
			out[out==" NA NA"] = "000000"
			for(c in 1:dim(out)[1]){
				cat(paste(out[c,], collapse=" "), file=file, append=T)
				cat("\n", file= file, append=T)
			}	
		}
	
	close(file)	



###############################################
## TMVP: Output juvs and adults from one focal pop: TopTwo on anemone >= 8cm vs. lower quartile)
## Uses Genepop format with two populations. Run through Formatomatic to get to TMVP

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
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
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"


# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
groups =c("Ad", "Juv")
pops = sort(unique(genowide$SiteNum))
for(j in 1:length(pops)){	
	a = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Ad")
	b = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Juv")
	print(paste("Pops ", pops[j], ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

# Write out focal pop
pop = 25
gens = c(0,1) # generation times of the groups
cols = grep("A1.", names(genowide))
loci = unlist(strsplit(names(genowide)[cols], "A1.", fixed=T))[seq(2,26, by=2)]
	
	# Write to file
	name = paste("Aclarkii_", Sys.Date(), "_TVMP", pop, ".gen", sep="")
	file = file(name)
	open(file, "w")

	# Header
	cat(paste("Title line: Aclarkii pops, ", pop, " Ads and Juvs", sep=""), file= file, sep=",", append=F)
	cat("\n", file= file, append=T)
	cat(paste(loci, collapse="\n"), file= file, append=T)
	cat("\n", file= file, append=T) 
	for(i in 1:length(groups)){
		k = which(genowide$SiteNum==pop & genowide$AdJuv == groups[i])
		if(length(k)>0){
			cat("Pop", file= file, sep="", append=T)
			cat("\n", file= file, append=T)
		
			# collapse locus alleles with ""
			out = paste("     ", rep(gens[i], length(k)), ",", sep="") # name each indiv according to population name
			out = cbind(out, paste(formatC(genowide[k,cols[1]], width=3, flag="0"), formatC(genowide[k,(cols[1]+1)], width=3, flag="0"), sep=""))
			for(c in cols[2:length(cols)]){
				out = cbind(out, paste(formatC(genowide[k,c], width=3, flag="0"), formatC(genowide[k,(c+1)], width=3, flag="0"), sep=""))
			}
			out[out==" NA NA"] = "000000"
			for(c in 1:dim(out)[1]){
				cat(paste(out[c,], collapse=" "), file=file, append=T)
				cat("\n", file= file, append=T)
			}	
		}
	}
	close(file)	
}



########################
# Output for create and FaMoz

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)

dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

	# add population info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

	# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
	widenames = names(genowide)
	i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
	genowide = genowide[,-i]
	
	# add Ad/juv info: Adults are top two on anemone if >=8cm, juvs are everything else
	genowide$AdJuv = 0
	i = genowide$TopTwo & genowide$Size >= 8 
	genowide$AdJuv[i] = 0 # Adults
	genowide$AdJuv[!i] = 1 # Juvs

	# remove pops 18 and 7?
	genowide = genowide[genowide$SiteNum != 7 & genowide$SiteNum != 18,]

	
	# Write to file
#	file = file(paste("Aclarkii_create_", Sys.Date(), ".csv", sep=""))
	file = file(paste("Aclarkii_create_", Sys.Date(), "no18or7.csv", sep=""))
	open(file, "w")
	

	i = grep("A1", names(genowide))
	j = names(genowide)[i[1]:(i[length(i)]+1)]
	j = paste(j, collapse=",")
	j = gsub(".", "_", j, fixed=T)
	outline = paste("Pop,Indiv,Cohort,",j,sep="") 
	cat(outline, file=file, append=T)
	cat("\n", file=file, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[,c(grep("SiteNum",widenames),grep("Sample",widenames),grep("AdJuv", widenames), grep("A[[:digit:]]",widenames))]
	write.table(outfile, file=file, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(file)	

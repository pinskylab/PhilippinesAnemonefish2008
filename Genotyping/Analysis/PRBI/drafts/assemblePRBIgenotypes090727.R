setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")

#### Read the genotypes ####

colstoread = c("Sample.Name", "Marker", "Allele.1", "Allele.2")


afg1 = read.table("AFG1-3_PRBI Genotypes 090703.txt", header = T, sep ="\t", na.string="?")
	afg1 = afg1[,colstoread]
afg9 = read.table("AFG9-10_PRBI Genotypes 090621.txt", header = T, sep ="\t", na.string="?")
	afg9 = afg9[,colstoread]

afg1_3.m.mlp = read.table("AFG1-3.M_PRBI_MLP Genotypes 090703.txt", header = T, sep ="\t", na.string="?")
	afg1_3.m.mlp = afg1_3.m.mlp[,colstoread]
afg4_5.m.mlp = read.table("AFG4-5.M_PRBI_MLP Genotypes 090709.txt", header = T, sep ="\t", na.string="?")
	afg4_5.m.mlp = afg4_5.m.mlp[,colstoread]
afg6_7.m.mlp = read.table("AFG6-7.M_PRBI_MLP Genotypes 090729.txt", header = T, sep ="\t", na.string="?")
	afg6_7.m.mlp = afg6_7.m.mlp[,colstoread]
afg8.m.mlp = read.table("AFG8.M_PRBI_MLP Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg8.m.mlp = afg8.m.mlp[,colstoread]
afg9.m.mlp = read.table("AFG9.M_PRBI_MLP Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg9.m.mlp = afg9.m.mlp[,colstoread]
afg45 = read.table("AFG45_PRBI Genotypes 090729.txt", header = T, sep ="\t", na.string="?")
	afg45 = afg45[,colstoread]
afg46_47.mlp = read.table("AFG46-47_PRBI_MLP Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg46_47.mlp = afg46_47.mlp[,colstoread]
afg48.mlp = read.table("AFG48_PRBI_MLP Genotypes 090724.txt", header = T, sep ="\t", na.string="?")
	afg48.mlp = afg48.mlp[,colstoread]
afg49.mlp = read.table("AFG49_PRBI_MLP Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg49.mlp = afg49.mlp[,colstoread]
afg50 = read.table("AFG50_PRBI_MLP Genotypes 090714.txt", header = T, sep ="\t", na.string="?")
	afg50 = afg50[,colstoread]
afg51_57.mlp = read.table("AFG51-57_PRBI_MLP Genotypes 090729.txt", header = T, sep ="\t", na.string="?")
	afg51_57.mlp = afg51_57.mlp[,colstoread]


afg1.m.mg = read.table("../MattGribble/AFG1.M_PRBI_MG Genotypes 090702.txt", header = T, sep ="\t", na.string="?")
	afg1.m.mg = afg1.m.mg[,colstoread]
afg2.m.mg = read.table("../MattGribble/AFG2.M_PRBI_MG Genotypes 090702.txt", header = T, sep ="\t", na.string="?")
	afg2.m.mg = afg2.m.mg[,colstoread]
afg3.m.mg = read.table("../MattGribble/AFG3.M_PRBI_MG Genotypes 090702.txt", header = T, sep ="\t", na.string="?")
	afg3.m.mg = afg3.m.mg[,colstoread]
afg4_5.m.mg = read.table("../MattGribble/AFG4and5.M_PRBI_MG Genotypes 090702.txt", header = T, sep ="\t", na.string="?")
	afg4_5.m.mg = afg4_5.m.mg[,colstoread]
afg6_7.m.mg = read.table("../MattGribble/AFG6-7.M_PRBI_MG Genotypes 090727.txt", header=T, sep="\t", na.string="?")
	afg6_7.m.mg = afg6_7.m.mg[,colstoread]
afg8.m.mg = read.table("../MattGribble/AFG8.M_PRBI_MG Genotypes 090724.txt", header=T, sep="\t", na.string="?")
	afg8.m.mg = afg8.m.mg[,colstoread]
afg9.m.mg = read.table("../MattGribble/AFG9.M_PRBI_MG Genotypes 090727.txt", header=T, sep="\t", na.string="?")
	afg9.m.mg = afg9.m.mg[,colstoread]
afg46.mg = read.table("../MattGribble/AFG46_PRBI_MG Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg46.mg = afg46.mg[,colstoread]
afg47.mg = read.table("../MattGribble/AFG47_PRBI_MG Genotypes 090724.txt", header = T, sep ="\t", na.string="?")
	afg47.mg = afg47.mg[,colstoread]
afg48.mg = read.table("../MattGribble/AFG48_PRBI_MG Genotypes 090724.txt", header = T, sep ="\t", na.string="?")
	afg48.mg = afg48.mg[,colstoread]
afg49.mg = read.table("../MattGribble/AFG49_PRBI_MG Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg49.mg = afg49.mg[,colstoread]

afg51.mg = read.table("../MattGribble/AFG51_PRBI_MG Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg51.mg = afg51.mg[,colstoread]
afg52.mg = read.table("../MattGribble/AFG52_PRBI_MG Genotypes 090723.txt", header = T, sep ="\t", na.string="?")
	afg52.mg = afg52.mg[,colstoread]
afg53.mg = read.table("../MattGribble/AFG53_PRBI_MG Genotypes 090724.txt", header = T, sep ="\t", na.string="?")
	afg53.mg = afg53.mg[,colstoread]
afg54.mg = read.table("../MattGribble/AFG54_PRBI_MG Genotypes 090727.txt", header = T, sep ="\t", na.string="?")
	afg54.mg = afg54.mg[,colstoread]
afg55.mg = read.table("../MattGribble/AFG55_PRBI_MG Genotypes 090724.txt", header = T, sep ="\t", na.string="?")
	afg55.mg = afg55.mg[,colstoread]
afg56.mg = read.table("../MattGribble/AFG56_PRBI_MG Genotypes 090729.txt", header = T, sep ="\t", na.string="?")
	afg56.mg = afg56.mg[,colstoread]
afg57.mg = read.table("../MattGribble/AFG57_PRBI_MG Genotypes 090723.txt", header = T, sep ="\t", na.string="?")
	afg57.mg = afg57.mg[,colstoread]


###########################################
#### Compare MLP vs. MG genotype calls ####

colnames = c("Sample", "Marker", "A1", "A2")
colnames(afg1.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg2.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg3.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg3.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg4_5.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg6_7.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg8.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg9.m.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg46.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg47.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg48.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg49.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg51.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg52.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg53.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg54.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg55.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg56.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")
colnames(afg57.mg) = paste(colnames, c("","",".MG", ".MG"), sep="")

colnames(afg1_3.m.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg4_5.m.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg6_7.m.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg8.m.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg9.m.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg46_47.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg48.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg49.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")
colnames(afg51_57.mlp) = paste(colnames, c("","",".MLP", ".MLP"), sep="")


compareAFG = function(mg, mlp){
	mg$SampleMarker = paste(mg$Sample, mg$Marker, sep="")
	mlp$SampleMarker = paste(mlp$Sample, mlp$Marker, sep="")

	# remove duplicates, preserving the one with genotypes
	i = order(mg$A1.MG, mg$A2.MG, na.last=TRUE)
	mg = mg[i,]
	mg = mg[!duplicated(mg$SampleMarker),]
	i = order(mlp$A1.MLP, mlp$A2.MLP, na.last=TRUE)
	mlp = mlp[i,]
	mlp = mlp[!duplicated(mlp$SampleMarker),]

	afg = merge(mg, subset(mlp, select=c(-Marker, -Sample)), by="SampleMarker")

	afg$A1 = NA
	afg$A2 = NA
	afg$error = NA
	for(i in 1:length(afg$SampleMarker)){
		a = identical(afg$A1.MG[i], afg$A1.MLP[i])
		b = identical(afg$A2.MG[i], afg$A2.MLP[i])
		if(a & b){
			afg$A1[i] = afg$A1.MLP[i]
			afg$A1[i] = afg$A1.MLP[i]	
			afg$error[i] = FALSE	
		} else {
			afg$error[i] = TRUE	
		}
	}

	print(paste("Errors: ", sum(afg$error), " (", round(sum(afg$error)/length(afg$error)*100, digits=2), "%)", sep=""))

	return(afg)
}


# AFG1-3.M 
	afg1_3.m.mg = rbind(afg1.m.mg, afg2.m.mg, afg3.m.mg)
	afg = compareAFG(afg1_3.m.mg, afg1_3.m.mlp)
	afg[afg$error,]

	# Choose the MLP genotype calls (were checked by eye!)
	afg1.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

	
## AFG4-5.M

	afg = compareAFG(afg4_5.m.mg, afg4_5.m.mlp)
	print(paste("Errors: ", sum(afg$error), " (", round(sum(afg$error)/length(afg$error)*100, digits=2), "%)", sep=""))
	print(paste("Errors: ", sum(afg$error[afg$Marker != "APY_2"]), " (", round(sum(afg$error[afg$Marker != "APY_2"])/length(afg$error[afg$Marker != "APY_2"])*100, digits=2), "%)", sep=""))

	afg[afg$error,]
	afg[afg$error & afg$Marker != "APY_2",]

	# choose MLP genotypes as consensus (only after checking by eye!)
	afg4.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]
	

## AFG6-7.M
	afg = compareAFG(afg6_7.m.mg, afg6_7.m.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]
	afg6.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AFG8.M
	afg = compareAFG(afg8.m.mg, afg8.m.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]
	afg8.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AF9.M
	afg = compareAFG(afg9.m.mg, afg9.m.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]
	afg9.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

# AFG46-47.M 
	afg46_47.mg = rbind(afg46.mg, afg47.mg)
	afg = compareAFG(afg46_47.mg, afg46_47.mlp)
	afg[afg$error,]

	# Choose the MLP genotype calls (were checked by eye!)
	afg46 = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AF48
	afg = compareAFG(afg48.mg, afg48.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]
	afg48 = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AF49
	afg = compareAFG(afg49.mg, afg49.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]
	afg49 = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

# AFG51-57.M 
	afg51_57.mg = rbind(afg51.mg, afg52.mg, afg53.mg, afg54.mg, afg55.mg, afg56.mg, afg57.mg)
	afg = compareAFG(afg51_57.mg, afg51_57.mlp)
	afg[afg$error,]

	# Choose the MLP genotype calls (were checked by eye!)
	afg51 = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

	
#### Assemble the genotypes ####

colnames = c("Sample", "Marker", "A1", "A2")
colnames(afg1) = paste(colnames, c("","",".AFG1", ".AFG1"), sep="")
colnames(afg9) = paste(colnames, c("","",".AFG9", ".AFG9"), sep="")
colnames(afg1.m) = paste(colnames, c("","",".AFG1.M", ".AFG1.M"), sep="")
colnames(afg4.m) = paste(colnames, c("","",".AFG4.M", ".AFG4.M"), sep="")
colnames(afg6.m) = paste(colnames, c("","",".AFG6.M", ".AFG6.M"), sep="")
colnames(afg8.m) = paste(colnames, c("","",".AFG8.M", ".AFG8.M"), sep="")
colnames(afg9.m) = paste(colnames, c("","",".AFG9.M", ".AFG9.M"), sep="")
colnames(afg45) = paste(colnames, c("","",".AFG45", ".AFG45"), sep="")
colnames(afg46) = paste(colnames, c("","",".AFG46", ".AFG46"), sep="")
colnames(afg48) = paste(colnames, c("","",".AFG48", ".AFG48"), sep="")
colnames(afg49) = paste(colnames, c("","",".AFG49", ".AFG49"), sep="")
colnames(afg50) = paste(colnames, c("","",".AFG50", ".AFG50"), sep="")
colnames(afg51) = paste(colnames, c("","",".AFG51", ".AFG51"), sep="")

		
# trim to PRBI samples and focal loci
afg1 = afg1[intersect(grep("PRBI",afg1$Sample), which(!is.na(afg1$A1))),] # remove extra loci from AFG1
afg9 = afg9[grep("PRBI",afg9$Sample),]
afg1.m = afg1.m[intersect(grep("PRBI",afg1.m$Sample), which(!is.na(afg1.m$A1))),] # remove extra genotypes from AFG1-3.M
afg4.m = afg4.m[intersect(grep("PRBI",afg4.m$Sample), which(!is.na(afg4.m$A1))),] # remove extra genotypes from AFG4 and 5.M
afg6.m = afg6.m[intersect(grep("PRBI",afg6.m$Sample), which(!is.na(afg6.m$A1))),] # remove extra genotypes from AFG6 and 7.M
afg8.m = afg8.m[intersect(grep("PRBI",afg8.m$Sample), which(!is.na(afg8.m$A1))),] # remove extra genotypes from AFG8.M
afg9.m = afg9.m[intersect(grep("PRBI",afg9.m$Sample), which(!is.na(afg9.m$A1))),] # remove extra genotypes from AFG9.M
afg45 = afg45[intersect(grep("PRBI",afg45$Sample), which(!is.na(afg45$A1))),] # remove extra genotypes from AFG45
afg46 = afg46[intersect(grep("PRBI",afg46$Sample), which(!is.na(afg46$A1))),] # remove extra genotypes from AFG46-47
afg48 = afg48[intersect(grep("PRBI",afg48$Sample), which(!is.na(afg48$A1))),] # remove extra genotypes from AFG48
afg49 = afg49[intersect(grep("PRBI",afg49$Sample), which(!is.na(afg49$A1))),] # remove extra genotypes from AFG49
afg50 = afg50[intersect(grep("PRBI",afg50$Sample), which(!is.na(afg50$A1))),] # remove extra genotypes from AFG50
afg51 = afg51[intersect(grep("PRBI",afg51$Sample), which(!is.na(afg51$A1))),] # remove extra genotypes from AFG51

#check for duplicated sample/loci combinations
which(duplicated(paste(afg1$Sample, afg1$Marker))) # yes
which(duplicated(paste(afg9$Sample, afg9$Marker)))
which(duplicated(paste(afg1.m$Sample, afg1.m$Marker)))
which(duplicated(paste(afg4.m$Sample, afg4.m$Marker)))
which(duplicated(paste(afg6.m$Sample, afg6.m$Marker)))
which(duplicated(paste(afg8.m$Sample, afg8.m$Marker)))
which(duplicated(paste(afg9.m$Sample, afg9.m$Marker)))
which(duplicated(paste(afg45$Sample, afg45$Marker)))
which(duplicated(paste(afg46$Sample, afg46$Marker)))
which(duplicated(paste(afg48$Sample, afg48$Marker)))
which(duplicated(paste(afg49$Sample, afg49$Marker)))
which(duplicated(paste(afg50$Sample, afg50$Marker)))
which(duplicated(paste(afg51$Sample, afg51$Marker)))

# do duplicates in files have identical genotypes?
# equalities should return TRUE if all genos are identical
sum(duplicated(paste(afg1$Sample, afg1$Marker))) == sum(duplicated(paste(afg1$Sample, afg1$Marker, afg1$A1, afg1$A2)))

# trim duplicates from file_ (AFG__) and file_ (AFG__)
afg1 = afg1[!duplicated(paste(afg1$Sample, afg1$Marker)),]

dim(afg1)
dim(afg9)
dim(afg1.m)
dim(afg4.m)
dim(afg6.m)
dim(afg8.m)
dim(afg9.m)
dim(afg45)
dim(afg46)
dim(afg48)
dim(afg49)
dim(afg50)
dim(afg51)
dim(afg52)
dim(afg53)
dim(afg54)
dim(afg55)
dim(afg56)
dim(afg57)


geno = merge(afg1, afg9, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg1.m, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg4.m, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg6.m, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg8.m, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg9.m, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg45, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg46, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg48, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg49, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg50, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg51, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)



# Fill in homozygotes
alleleones = grep("A1.", names(geno), fixed=T)
alleletwos = grep("A2.", names(geno), fixed=T)
for(i in 1:length(alleleones)){
	j <- !is.na(geno[,alleleones[i]]) & is.na(geno[,alleletwos[i]]) # find homozygotes
	geno[j,alleletwos[i]] <- geno[j,alleleones[i]]
}

# Make sure we list each individual for each locus in geno, remove other loci
remove = c("APY_2", "APY_44", "APY_45", "APY_65")
for(i in 1:length(remove)){
	j = which(geno$Marker==remove[i])
	geno = geno[-j,]	
}
markers = unique(geno$Marker)
markers = markers[!is.na(markers)]
geno= geno[geno$Sample!="PRBIContam",] # remove a contaminated well
indivs = unique(geno$Sample)
for(i in 1:length(markers)){	
	theseindivs = unique(geno$Sample[geno$Marker==markers[i]])
	if(length(theseindivs) < length(indivs)){
		new = setdiff(indivs, theseindivs)
		add = data.frame(Sample = new, Marker = markers[i])
		geno = merge(geno, add, all.x=T, all.y=T)
		dim(geno)
	}	
	if(length(theseindivs) > length(indivs)){
		print(paste("Too many indivs for marker ", marker[i]))
	}
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

# Find errors
err=table(geno$error, geno$Marker)
if(dim(err)[1]==1){
	err = rbind(err, rep(0, dim(err)[2]))	
	rownames(err) = c("FALSE", "TRUE")
}
err = rbind(err, (err[1,]+err[2,]))
err = rbind(err,round((err[2,]/err[3,])*100,1))
rownames(err) = c("No Errors", "Errors", "# Checked", "Error %")
t(err)

# Write out errors
i = which(geno$error==T)
write.csv(geno[i,], paste("Genotyping_errors_PRBI", Sys.Date(), ".csv", sep=""))

# Select a subset to redo for errorchecking (ran all on 7/13/09)
#redo = data.frame(Sample=character(0), Marker=character(0))
#redomarks = markers[markers!="AC1359" & markers!="AC1578" & markers!="APR_Cf29" & markers!="NNG_012"] # was 7/13/09
#for(i in 1:length(redomarks)){
#	a = err[3, which(colnames(err)==redomarks[i])]
#	if(a < 20){
#		redo = rbind(redo, data.frame(Sample=sample(unique(geno$Sample), 20), Marker=redomarks[i]))
#	} 	
#}
#dim(redo)
#write.csv(redo, paste("RedoPRBI_", Sys.Date(), ".csv", sep=""), row.names=FALSE)



# Calc % complete
i = table(!is.na(geno$A1consens), geno$Marker)
i = rbind(i,round((i[2,]/(i[1,]+i[2,]))*100,1))
rownames(i) = c("Missing", "Complete", "%")
t(i)
colSums(i[1:2,]) # number of indivs at each locus

# Which indivs are still missing at each locus? Write that out if desired
out = data.frame(Marker = character(0), Sample = character(0))
#missingmarkers = c("ACH_A11", "ACH_A4", "ACH_B9", "NNG_028") # was for 7/13/09
missingmarkers = c("ACH_A7", "ACH_C1", "ACH_D1", "ACH_A3", "ACH_A8", "NNG_004", "NNG_007")
for(i in 1:length(missingmarkers)){
	a = geno[geno$Marker==missingmarkers[i] & is.na(geno$A1consens), c("Marker", "Sample")]
	out = rbind(out, a)
}
out$reason = "missing"
dim(out)
#write.csv(out, paste("MissingPRBI_", Sys.Date(), ".csv", sep=""), row.names=FALSE)



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
for(i in 1:length(loci)){
	j=sort(unique(c(geno$A1consens[geno$Marker == loci[i]], geno$A2consens[geno$Marker== loci[i]])))
	cat(paste(loci[i],paste(j, collapse=" ")))
	cat("\n")
}


###################################################
# Write out whole table
write.csv(geno, paste("PRBI_genotypes_", Sys.Date(), ".csv", sep=""))
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
	
## REMOVE APCL285, 286 because identical
i = genowide$Sample == "APCL285" | genowide$Sample == "APCL286"
genowide = genowide[-i,]	
	
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
geno = read.csv("PRBI_genotypes_2009-07-28.csv", row.names=1)

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
	
	# add Ad/juv info: Adults are top two on anemone if >=8cm, juvs are everything else
	cutoff = 8
	genowide$AdJuv = 0
	i = genowide$TopTwo & genowide$Size > cutoff
	genowide$AdJuv[i] = 0 # Adults
	genowide$AdJuv[!i] = 1 # Juvs

	# remove pops 18 and 7?
	#genowide = genowide[genowide$SiteNum != 7 & genowide$SiteNum != 18,]

	
	# Write to file
	file = file(paste("PRBI_create_cut", cutoff, "_", Sys.Date(), ".csv", sep=""))
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

	filename = paste("PRBI_genowide_", Sys.Date(), ".csv", sep="")
	# write.csv(genowide, filename)
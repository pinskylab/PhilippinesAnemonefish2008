setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")

#### Read the genotypes ####

colstoread = c("Sample.Name", "Marker", "Allele.1", "Allele.2")

# only malin
afg1 = read.table("AFG1-3_PRBI Genotypes 090729.txt", header = T, sep ="\t", na.string="?")
	afg1 = afg1[,colstoread]
afg9 = read.table("AFG9-10_PRBI Genotypes 090621.txt", header = T, sep ="\t", na.string="?")
	afg9 = afg9[,colstoread]
afg15 = read.table("AFG15_PRBI Genotypes 090623.txt", header=T, sep="\t", na.string="?")
	afg15 = afg15[,colstoread]
afg58 = read.table("AFG58-59_PRBI_MLP Genotypes 090830.txt", header = T, sep ="\t", na.string="?")
	afg58 = afg58[,colstoread]
afg60 = read.table("AFG60-61_PRBI_MLP Genotypes 090922.txt", header = T, sep ="\t", na.string="?")
	afg60 = afg60[,colstoread]
b2 = read.table("B2-4 to C2-6_PRBI_MLP Genotypes 091111.txt", header=T, sep="\t", na.string="?")
	b2 = b2[,colstoread]

# malin to match against matt's calls
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

# matt's calls to match against malin
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

	# Choose the MLP genotype calls (were checked by eye!)
	afg6.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AFG8.M
	afg = compareAFG(afg8.m.mg, afg8.m.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]

	# Choose the MLP genotype calls (were checked by eye!)
	afg8.m = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AF9.M
	afg = compareAFG(afg9.m.mg, afg9.m.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]

	# Choose the MLP genotype calls (were checked by eye!)
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

	# Choose the MLP genotype calls (were checked by eye!)
	afg48 = afg[,c("Sample", "Marker", "A1.MLP", "A2.MLP")]

## AF49
	afg = compareAFG(afg49.mg, afg49.mlp)
	afg[afg$error,]
	i = afg$error
	afg[i,][order(afg$Marker[i], afg$Sample[i]),]

	# Choose the MLP genotype calls (were checked by eye!)
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
colnames(afg15) = paste(colnames, c("","",".AFG15", ".AFG15"), sep="")
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
colnames(afg58) = paste(colnames, c("","",".AFG58", ".AFG58"), sep="")
colnames(afg60) = paste(colnames, c("","",".AFG60", ".AFG60"), sep="")
colnames(b2) = paste(colnames, c("","",".B2", ".B2"), sep="")
		
# trim to PRBI samples and focal loci
afg1 = afg1[intersect(grep("PRBI",afg1$Sample), which(!is.na(afg1$A1))),] # remove extra loci from AFG1
afg9 = afg9[grep("PRBI",afg9$Sample),]
afg15 = afg15[intersect(grep("PRBI",afg15$Sample), which(!is.na(afg15$A1))),] # remove extra loci from AFG15
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
afg58 = afg58[intersect(grep("PRBI",afg58$Sample), which(!is.na(afg58$A1))),] # remove extra genotypes from AFG58
afg60 = afg60[intersect(grep("PRBI",afg60$Sample), which(!is.na(afg60$A1))),] # remove extra genotypes from AFG60
b2 = b2[intersect(grep("PRBI",b2$Sample), which(!is.na(b2$A1))),] # remove extra genotypes from B2-4 to C2-6

#check for duplicated sample/loci combinations
which(duplicated(paste(afg1$Sample, afg1$Marker))) # yes
which(duplicated(paste(afg9$Sample, afg9$Marker)))
which(duplicated(paste(afg15$Sample, afg15$Marker)))
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
which(duplicated(paste(afg58$Sample, afg58$Marker)))
which(duplicated(paste(afg60$Sample, afg60$Marker)))
which(duplicated(paste(b2$Sample, b2$Marker)))

# do duplicates in files have identical genotypes?
# equalities should return TRUE if all genos are identical
sum(duplicated(paste(afg1$Sample, afg1$Marker))) == sum(duplicated(paste(afg1$Sample, afg1$Marker, afg1$A1, afg1$A2)))

# trim duplicates (only if duplicates have identical genotypes!)
afg1 = afg1[!duplicated(paste(afg1$Sample, afg1$Marker)),]

dim(afg1)
dim(afg9)
dim(afg15)
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
dim(afg58)
dim(afg60)
dim(b2)

geno = merge(afg1, afg9, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg15, all.x=T, all.y = T, by = c("Sample", "Marker"))
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
geno = merge(geno, afg58, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, afg60, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)
geno = merge(geno, b2, all.x=T, all.y = T, by = c("Sample", "Marker"))
dim(geno)

# Remove AFG58-59 genotypes for NNG_028 or NNG_004 since these were Fermentas and don't count as a verification
i = geno$Marker == "NNG_028" | geno$Marker=="NNG_004"
geno$A1.AFG58[i] = NA
geno$A2.AFG58[i] = NA

# Fill in homozygotes
alleleones = grep("A1.", names(geno), fixed=T)
alleletwos = grep("A2.", names(geno), fixed=T)
for(i in 1:length(alleleones)){
	j <- !is.na(geno[,alleleones[i]]) & is.na(geno[,alleletwos[i]]) # find homozygotes
	geno[j,alleletwos[i]] <- geno[j,alleleones[i]]
}

# Remove loci we don't want
remove = c("APY_2", "APY_44", "APY_45", "APR_Cf8", "ACH_D8", "APY_65")
dim(geno)
for(i in 1:length(remove)){
	j = which(geno$Marker==remove[i])
	geno = geno[-j,]	
}
dim(geno)


# Make sure we list each individual for each locus in geno
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


# Order geno by Marker and Sample
i = order(geno$Marker, geno$Sample) # order by marker then sample
geno = geno[i,]
geno$Marker = as.character(geno$Marker)


# Find errors
err=table(geno$error, geno$Marker)
if(dim(err)[1]==1){
	err = rbind(err, rep(0, dim(err)[2]))	
	rownames(err) = c("FALSE", "TRUE")
}
err = rbind(err, (err[1,]+err[2,]))
err = rbind(err,round((err[2,]/err[3,])*100,1))
rownames(err) = c("No Errors", "Errors", "# Checked", "Error %")
err = t(err[,order(colnames(err))])
err

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
i = t(i[,order(colnames(i))])
i
colSums(i[1:2,]) # number of indivs at each locus

# Which indivs are still missing at each locus? Write that out if desired
out = data.frame(Marker = character(0), Sample = character(0))
#missingmarkers = c("ACH_A11", "ACH_A4", "ACH_B9", "NNG_028") # was for 7/13/09
missingmarkers = c("ACH_A7", "ACH_C1", "ACH_D1", "ACH_A3", "ACH_A8", "NNG_004", "NNG_007")
missingmarkers = c("AC1359", "AC1578", "ACH_A11", "ACH_A3", "ACH_A4", "ACH_A7",  "ACH_A8", "ACH_B9", "ACH_C1", "ACH_D1", "APR_Cf29", "APR_Cf39", "NNG_004", "NNG_007", "NNG_012", "NNG_028")
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
j = matrix(data= seq(0,13), nrow=3, ncol=14, byrow=T)
j[2,] = hist(i[1,], plot=F, breaks = seq(-0.5, 13.5, by=1))$counts
j[3,] = round(j[2,]/160*100)
rownames(j) = c("Num Missing", "Num Samples", "% Samples")
print(j)

# Remove indivs where % missing is very high
i = table(!is.na(geno$A1consens), geno$Sample)
hist(i[1,], plot=F, breaks = seq(-0.5, max(i[1,])+2, by=1))$counts
j = which(i[1,]>5)
print(k<-colnames(i)[j])
j=numeric(0)
for(i in k){
	j = c(j, which(geno$Sample==i))
}
dim(geno)
unique(geno$Sample[j])
geno = geno[-j,]
dim(geno)

### Fix genotyping errors where calls were fixable

# Check for allele dropout and choose genotype with more alleles as long as at least one allele matches
err = which(geno$error)
alleleones = grep("A1.", names(geno), fixed=T)
alleletwos = grep("A2.", names(geno), fixed=T)
for(i in err){
	checkone = alleleones[!is.na(geno[i,alleleones])] # get allele cols that aren't na
	checktwo = alleletwos[!is.na(geno[i,alleletwos])] # get allele cols that aren't na

	if(length(checkone) != length(checktwo)){ print(paste("A1 and A2 not equal length", i))}

	het = rep(NA, length(checkone))
	for(j in 1:length(checkone)){ # for each genotype
		het[j] = geno[i,checkone[j]] != geno[i,checktwo[j]] # which genotypes are heterozygotes?
	}
	if(sum(het)==1){ # if one het and rest hom, check that homs match and match the het, then fill in consensus 
		homsmatch = FALSE # flag as to whether or not the homs match themselves
		hommatchhets = FALSE # flag as to whether or not the homs match the het
		homs = which(!het)
		if(length(homs)>1){
			if(all(geno[i,checkone[homs]] == geno[i,checkone[homs][1]])) homsmatch = TRUE # all homs are the same
		} else {
			homsmatch = TRUE # only one hom, so of course it matches itself
		}
		if(homsmatch){ # if all the homs are the same
			m1 = geno[i,checkone[homs][1]] == geno[i,checkone[het][1]] # hom matches first allele of het?
			m2 = geno[i,checktwo[homs][1]] == geno[i,checktwo[het][1]] # hom matches second allele of het?
			if(m1 | m2){ # if hom allele matches one of the het alleles
				hommatchhets = TRUE
			}
		}
		if(homsmatch & hommatchhets){ # if both checks are true, fill in consensus genotypes
			geno$A1consens[i] = geno[i,checkone[het][1]]
			geno$A2consens[i] = geno[i,checktwo[het][1]]
			print(paste("Fixed", geno$Sample[i], geno$Marker[i]))
		} else {
			print(paste("FAILED: Didn't match homs to hets for", geno$Sample[i], geno$Marker[i]))	
		}
	}
	if(sum(het)>1){ # if > one het, check for equality of the hets and homs first
		homsmatch = FALSE # flag as to whether or not the homs match themselves
		hommatchhets = FALSE # flag as to whether or not the homs match the het
		homs = which(!het)
		if(length(homs)>1){
			if(all(geno[i,checkone[homs]] == geno[i,checkone[homs][1]])) homsmatch = TRUE # all homs are the same
		} else {
			homsmatch = TRUE # only one hom, so of course it matches itself
		}
		if(homsmatch){ # if all the homs are the same
			m1 = geno[i,checkone[homs][1]] == geno[i,checkone[het][1]] # hom matches first allele of het?
			m2 = geno[i,checktwo[homs][1]] == geno[i,checktwo[het][1]] # hom matches second allele of het?
			if(m1 | m2){ # if hom allele matches one of the het alleles
				hommatchhets = TRUE
			}
		}
		matchone = all(geno[i,checkone[het]] == geno[i,checkone[which(het)[1]]]) # check for equality of hets
		matchtwo = all(geno[i,checktwo[het]] == geno[i,checktwo[which(het)[1]]]) # check for equality of hets
		if(matchone & matchtwo & homsmatch & hommatchhets){
			geno$A1consens[i] = geno[i,checkone[het][1]]
			geno$A2consens[i] = geno[i,checktwo[het][1]]
			print(paste("Fixed", geno$Sample[i], geno$Marker[i]))
		} else {
			print(paste("Hets and homs don't match for", geno$Sample[i], geno$Marker[i]))
		}
	}	
	
	if((i %% 500) == 0){ print(i)}
		
}

# Check for homozygotes
geno$hom = geno$A1consens == geno$A2consens
#geno$hom[is.na(geno$hom)] = FALSE
hom = table(geno$Marker, geno$hom)
hom = cbind(hom, hom[,1]+hom[,2], round((hom[,2]/(hom[,1]+hom[,2]))*100,1))
hom = hom[order(rownames(hom)),]
colnames(hom) = c("Hets", "Homs", "Total", "%Hom")
hom

# Find unverified homozygotes for ACH_B9, NNG_004, and NNG_007
# Acceptable locus/Taq combinations:
#	ACH_B9 / Fermentas: AFG53 (in 51-57), 58, 59 (in 58-59), 60, 61 (in 60-61)
#	NNG_007 / Fermentas: AFG54 & 57 (in 51-57), 58-59, 60-61
#	NNG_004 / TypeIt: AFG49, 9.M, 60-61
#	NNG_028 / TypeIt: AFG6-7.M, 60-61
i = geno$hom & geno$Marker == "ACH_B9" & is.na(geno$A1.AFG51) & is.na(geno$A1.AFG58) & is.na(geno$A1.AFG60) & is.na(geno$A1.B2)
i2 = geno$hom & geno$Marker == "NNG_007" & is.na(geno$A1.AFG51) & is.na(geno$A1.AFG58) & is.na(geno$A1.AFG60) & is.na(geno$A1.B2)
i3 = geno$hom & geno$Marker == "NNG_004" & is.na(geno$A1.AFG49) & is.na(geno$A1.AFG9.M) & is.na(geno$A1.AFG60)
i4 = geno$hom & geno$Marker == "NNG_028" & is.na(geno$A1.AFG6.M) & is.na(geno$A1.AFG60)

j = i|i2|i3|i4
redos = data.frame(Sample=geno$Sample[j], ACH_B9=i[j], NNG_007=i2[j], NNG_004=i3[j], NNG_028=i4[j])
redos


# Calc % complete now that we've filled in the ones that we could fix
i = table(!is.na(geno$A1consens), geno$Marker)
i = rbind(i,round((i[2,]/(i[1,]+i[2,]))*100,1))
rownames(i) = c("Missing", "Complete", "%")
i = t(i[,order(colnames(i))])
i


# Check loci for regular alleles (no alleles 1bp apart, should match bins in GeneMapper)
loci = unique(geno$Marker)
for(i in 1:length(loci)){
	j=sort(unique(c(geno$A1consens[geno$Marker == loci[i]], geno$A2consens[geno$Marker== loci[i]])))
	cat(paste(loci[i],paste(j, collapse=" ")))
	if(length(j)>1){
		for(k in 2:length(j)){
			if(j[k] - j[k-1] == 1) print("ERROR")		}
	}
	cat("\n")
}

# Compare against the previous genotypes: what changed?
old = read.csv("PRBI_genotypes_2009-09-23.csv", row.names=1)
length(unique(geno$Sample))
length(unique(old$Sample))
all(unique(as.character(geno$Sample)) == unique(as.character(old$Sample))) # do all sample names match?
match = rep(0, dim(geno)[1])
for(i in 1:dim(geno)[1]){
	m1 = FALSE # alleles 1 match?
	m2 = FALSE # alleles 2 match?
	if(is.na(geno$A1consens[i]) & is.na(old$A1consens[i])){m1 = TRUE} # if both are NA
	if(is.na(geno$A2consens[i]) & is.na(old$A2consens[i])) {m2 = TRUE}

	if(xor(is.na(geno$A1consens[i]), is.na(old$A1consens[i]))){m1 = FALSE} # if one is NA
	if(xor(is.na(geno$A2consens[i]), is.na(old$A2consens[i]))){m1 = FALSE} # if one is NA

	if(!is.na(geno$A1consens[i]) & !is.na(old$A1consens[i])){ # if both are not NA
		m1 = geno$A1consens[i] == old$A1consens[i]
	}
	if(!is.na(geno$A2consens[i]) & !is.na(old$A2consens[i])){
		m2 = geno$A2consens[i] == old$A2consens[i]
	}
	match[i] = m1 & m2
}
i = which(match==0)
cbind(geno[i,c("Marker", "A1consens", "A2consens")], old[i,c("A1consens", "A2consens")])




###################################################
# Write out whole table
write.csv(geno, paste("PRBI_genotypes_", Sys.Date(), ".csv", sep=""))
###################################################

###################################################
# Read in whole table (geno)

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("PRBI_genotypes_2009-11-12.csv", row.names=1)
###################################################

#######################
# Make genowide

dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos", "hom")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,lat,long, Size)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

	# add UTM coords
	library(PBSmapping)
	lldata = genowide[,c("lat", "long")]
	names(lldata) = c("Y", "X")
	attr(lldata, "projection") <- "LL"
	utmdata = convUL(lldata, km=FALSE)
	names(utmdata) = c("y_utm", "x_utm")
	genowide = cbind(genowide, utmdata)

	# sort by region and sample ID
	i = order(genowide$Region, genowide$SiteNum, genowide$Sample)
	genowide = genowide[i,]
	
write.csv(genowide, paste("PRBI_genowide_", Sys.Date(), ".csv", sep=""))




# Sample sizes by site
table(genowide$SiteNum)


# GO TO formatPRBIgenotypes090831.R FOR OUTPUTS TO DIFFERENT PROGRAM FORMATS

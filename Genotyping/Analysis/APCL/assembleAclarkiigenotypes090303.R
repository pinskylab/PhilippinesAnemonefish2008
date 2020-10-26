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

# trim duplicates from file4 (AFG16)
file4 = file4[!duplicated(paste(file4$Sample, file4$Marker)),]

# update APY_45 allele calls in AFG17,18,20,21 to bp
file5$A1.AFG17[file5$Marker=="APY_45"] = (file5$A1.AFG17[file5$Marker=="APY_45"]-10)*2+210
file5$A2.AFG17[file5$Marker=="APY_45"] = (file5$A2.AFG17[file5$Marker=="APY_45"]-10)*2+210
	
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


# Write out whole table
write.csv(geno, paste("Aclarkii_genotypes_", Sys.Date(), ".csv", sep=""))
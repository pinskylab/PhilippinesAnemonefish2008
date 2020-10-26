setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Extractions")

toextract = read.csv("To Extract 081208.csv")
afe1 = as.matrix(read.csv("AFE1.csv", row.names=1))
	dim(afe1) = c(dim(afe1)[1]*dim(afe1)[2],1)

sum(duplicated(toextract$ID))
# == 0. Good!

# Remove AFE1 plate from our samples
temp = toextract[-match(afe1, toextract$ID),]


# Sample AFE2 so it complements AFE1 (17 APCL, 79 PRBI)
plate2 = sample(temp$ID[temp$Spp=="APCL"], 17, replace=F)
plate2 = c(plate2, sample(temp$ID[temp$Spp=="PRBI"], 79, replace=F))

temp2 = temp[-match(plate2, temp$ID),]
dim(plate2) = c(8,12)


# Sample AFE3 so that it's all APCL
plate3 = sample(temp2$ID[temp2$Spp=="APCL"], 96, replace=F)

temp3 = temp2[-match(plate3, temp2$ID),]
dim(plate3) = c(8,12)


# Sample AFE4 so that it's all APCL
plate4 = sample(temp3$ID[temp3$Spp=="APCL"], 96, replace=F)

remaining = temp3[-match(plate4, temp3$ID), ]
dim(plate4) = c(8,12)

# shouldn't be any matches
match(afe1, plate2)
match(afe1, plate3)
match(afe1, plate4)
match(afe1, remaining)
match(plate2, plate3)
match(plate2, plate4)
match(plate2, remaining)
match(plate3, plate4)
match(plate3, remaining)
match(plate4, remaining)

# should be the same
sum(length(afe1), length(plate2), length(plate3), length(plate4), length(remaining$ID))
length(toextract$ID)


# Write out results
write.csv(plate2, "AFE2.csv")
write.csv(plate3, "AFE3.csv")
write.csv(plate4, "AFE4.csv")
write.csv(remaining, "Remaining.csv")
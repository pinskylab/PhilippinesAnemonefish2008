setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Extractions")

toextract = read.csv("To Extract 081208.csv")
#write.csv(toextract, "To Extract 081208.csv")

sum(duplicated(toextract$ID))
# == 0. Good!

plate1 = sample(toextract$ID, 96, replace=F)

temp = toextract$ID[-match(plate1, toextract$ID)]

plate2 = sample(temp, 96, replace=F)

temp = temp[-match(plate2, temp)]

plate3 = sample(temp, 96, replace=F)

temp = temp[-match(plate3, temp)]

plate4 = sample(temp, 96, replace=F)

remaining = temp[-match(plate4, temp)]

# shouldn't be any matches
match(plate1, plate2)
match(plate1, plate3)
match(plate1, plate4)
match(plate1, remaining)
match(plate2, plate3)
match(plate2, plate4)
match(plate2, remaining)
match(plate3, plate4)
match(plate3, remaining)
match(plate4, remaining)

# should be the same
sum(length(plate1), length(plate2), length(plate3), length(plate4), length(remaining))
length(toextract$ID)


# Write out results
dim(plate1) = c(8,12)
dim(plate2) = c(8,12)
dim(plate3) = c(8,12)
dim(plate4) = c(8,12)

write.csv(plate1, "AFE1.csv")
write.csv(plate2, "AFE2.csv")
write.csv(plate3, "AFE3.csv")
write.csv(plate4, "AFE4.csv")
write.csv(remaining, "Remaining.csv")
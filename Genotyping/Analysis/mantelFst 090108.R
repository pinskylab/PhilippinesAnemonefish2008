setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008")

fst = read.csv("Genotyping/AFG5,6,7,8 fst matrix AC1359 090108.csv", row.names=1)
geo = read.csv("Genotyping/AFG5,6,7,8 dist matrix 090108.csv", row.names=1)

# linearize fst and set negs to zeros
fst[fst < 0 & !is.na(fst)] = 0
fstlin = fst/(1-fst)

cebu_fst = fstlin[2:11,2:11]
cebu_geo = geo[2:11,2:11]

leyte_fst = fstlin[12:20,12:20]
leyte_geo = geo[12:20, 12:20]

library(vegan)



mantel(cebu_geo, cebu_fst, method="pearson")

mantel(leyte_geo, leyte_fst, method="pearson")

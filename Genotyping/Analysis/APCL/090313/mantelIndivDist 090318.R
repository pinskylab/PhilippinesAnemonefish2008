setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008/Analysis/090313")

collect = read.csv("../../Collections2008-11-21.csv", header=T, row.names=1)
surveys = read.csv("../../GPSSurveys2009-01-08.csv", header=T)

geo = as.matrix(read.csv("Aclarkii_2009-03-13 indivgeo.csv", skip=2, header=T, row.names=379))
rouss = as.matrix(read.csv("Aclarkii_2009-03-13 rousseta.csv", skip=2, header=T, row.names=1))
rouss[upper.tri(rouss)]=NA

collect$sample = paste(collect$Spp, collect$ID, sep="")
meta = data.frame(sample = row.names(rouss))
meta$sample = as.character(meta$sample)
meta = merge(meta, subset(collect, select=c(sample, SurveyNum)), all.x=T, sort=F)
meta = merge(meta, subset(surveys, select=c(SurveyNum, Region)), all.x=T, sort=F)

identical(row.names(geo), row.names(rouss)) # must be TRUE!
identical(meta$sample, row.names(rouss)) # must be TRUE!

len = length(geo)

geolin = geo[lower.tri(geo)]
rousslin = rouss[lower.tri(rouss)]

length(geolin)
length(rousslin)

plot(geolin, rousslin)

l = lm(rousslin ~geolin)
summary(l)

library(vegan)
m = mantel(geo, rouss)


# Subset by Cebu or Leyte

i = (meta$Region == "Cebu")

c = rouss[i,i]
cg = geo[i,i]

plot(cg[lower.tri(cg)], c[lower.tri(c)])

i = (meta$Region == "Leyte")
l = rouss[i,i]
lg = geo[i,i]

plot(lg[lower.tri(lg)], l[lower.tri(l)])
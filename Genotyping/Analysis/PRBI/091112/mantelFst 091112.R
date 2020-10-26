library(vegan)

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/PRBI/091112/")

## Read in data

fst = read.csv("Genepop/PRBI_2009_11_12_fsts.csv", row.names=1) # Genepop with all pops
geo = read.csv("PRBI_2009-11-12_geo.csv", row.names=1) 

## The APCL data for comparison
apclfst = read.csv("../../APCL/090507/genepop/Aclarkii_091028_all_1pop_fst.csv", header=TRUE, row.names=1)
apclgeo = read.csv("../../APCL/090507/Aclarkii_2009-11-12 geo.csv", header=TRUE, row.names=1)


## Set up PRBI data
# linearize fst, remove diagonals
geo[upper.tri(geo, diag=T)] = NA
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)

# Drop Pop 15 (n=1)
i = grep("X15", names(fstlin))
if(length(i)>0) fstlinno15 = fstlin[-i, -i]
i = grep("X15", names(geo))
if(length(i)>0) geono15 = geo[-i, -i]

# Trim for comparison against APCL
fstlintrim = fstlinno15[2:11, 2:11]
geotrim = geono15[2:11, 2:11]


# Set up Cebu
cebu_fst = as.matrix(fstlin[3:10, 3:10])
cebu_geo = as.matrix(geo[3:10, 3:10])
cebu_fstno15 = as.matrix(fstlinno15[3:9, 3:9])
cebu_geono15 = as.matrix(geono15[3:9,3:9])

# Set up Leyte
leyte_fst = as.matrix(fstlinno15[10:11,10:11])
leyte_geo = as.matrix(geono15[10:11, 10:11])


## Set up APCL data
apclgeo[upper.tri(apclgeo, diag=T)] = NA
apclfst[upper.tri(apclfst, diag=T)] = NA
apclfstlin = apclfst/(1-apclfst)


# Drop Pops 15-18, 20, and 23-27 (to compare against prbi)
i = grep("X15|X16|X17|X18|X20|X23|X24|X25|X27", names(apclfst))
if(length(i)>0) apclfstlintrim = apclfstlin[-i, -i]

i = grep("X15|X16|X17|X18|X20|X23|X24|X25|X27", names(apclgeo))
if(length(i)>0) apclgeotrim = apclgeo[-i, -i]



apclcebu_fst = as.matrix(apclfstlin[2:11,2:11])
apclcebu_geo = as.matrix(apclgeo[2:11,2:11])
apclcebu_fsttrim = as.matrix(apclfstlintrim[2:8, 2:8]) # only the cebu pops that match PRBI
apclcebu_geotrim = as.matrix(apclgeotrim[2:8, 2:8])

apclleyte_fst = as.matrix(apclfstlin[12:19,12:19])
apclleyte_geo = as.matrix(apclgeo[12:19, 12:19])


#### Analysis

# Cebu and Leyte together
# Use pop 15
par(cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
xlims = c(0,250)
ylim = c(min(leyte_fst, cebu_fst, na.rm=T), max(leyte_fst, cebu_fst, na.rm=T))
col1 = "black" # Cebu
col2 = "grey60" # Leyte
cols = c(rep(col1, sum(!is.na(cebu_fstno15))), rep(col2, sum(!is.na(leyte_fst))))

plot(as.dist(cebu_geo), as.dist(cebu_fst), xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI Isolation-by-distance (w/ 15)", ylim=ylim, xlim = xlims, pch = 20, lwd=2, col=col1)
points(as.dist(leyte_geo), as.dist(leyte_fst), pch=20, lwd=2, col=col2)
legend("bottomright", legend=c("Cebu", "Leyte"), pch=20, col=c(col1, col2))
x = as.dist(cebu_geo)
y = as.dist(cebu_fst)
l=(lm(y ~ x))
lines(as.dist(cebu_geo)[!is.na(as.dist(cebu_geo))], l$fitted.values, col="blue")



# Cebu and Leyte together
# Don't use pop 15
par(cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
xlims = c(0,250)
ylim = c(min(leyte_fst, cebu_fstno15, na.rm=T), max(leyte_fst, cebu_fstno15, na.rm=T))
col1 = "black" # Cebu
col2 = "grey60" # Leyte
cols = c(rep(col1, sum(!is.na(cebu_fstno15))), rep(col2, sum(!is.na(leyte_fst))))

plot(as.dist(cebu_geono15), as.dist(cebu_fstno15), xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI Isolation-by-distance (no 15)", ylim=ylim, xlim = xlims, pch = 20, lwd=2, col=col1)
points(as.dist(leyte_geo), as.dist(leyte_fst), pch=20, lwd=2, col=col2)
legend("bottomright", legend=c("Cebu", "Leyte"), pch=20, col=c(col1, col2))
x = as.dist(cebu_geono15)
y = as.dist(cebu_fstno15)
l=(lm(y ~ x))
lines(as.dist(cebu_geono15)[!is.na(as.dist(cebu_geono15))], l$fitted.values, col="blue")



# Cebu and Leyte together
# Don't use pop 14 or 15
par(cex=1.4, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
xlims = c(0,250)
ylim = c(min(leyte_fst, cebu_fstno15[1:6,1:6], na.rm=T), max(leyte_fst, cebu_fstno15[1:6,1:6], na.rm=T))
col1 = "black" # Cebu
col2 = "grey60" # Leyte
cols = c(rep(col1, sum(!is.na(cebu_fstno15))), rep(col2, sum(!is.na(leyte_fst))))

plot(as.dist(cebu_geono15[1:6,1:6]), as.dist(cebu_fstno15[1:6,1:6]), xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI Isolation-by-distance (no 14 or 15)", ylim=ylim, xlim = xlims, pch = 20, lwd=2, col=col1)
points(as.dist(leyte_geo), as.dist(leyte_fst), pch=20, lwd=2, col=col2)
legend("bottomright", legend=c("Cebu", "Leyte"), pch=20, col=c(col1, col2))
x = as.dist(cebu_geono15[1:6,1:6])
y = as.dist(cebu_fstno15[1:6,1:6])
l=(lm(y ~ x))
lines(as.dist(cebu_geono15)[!is.na(as.dist(cebu_geono15))], l$fitted.values, col="blue")






# PRBI vs. APCL: Cebu slopes
quartz(width=9, height=5)
par(mfrow=c(1,2), cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
pchs = 20
lwds = 2
xlims = c(0,250)
ylim = c(min(c(apclcebu_fst, cebu_fst), na.rm=T), max(c(apclcebu_fst, cebu_fst), na.rm=T))
plot(cebu_geo[1:64], cebu_fst[1:64], pch=pchs, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI Cebu", ylim=ylim, xlim = xlims)
l = lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)])
lines(cebu_geo[1:64][!is.na(cebu_geo[1:64])], l$fitted.values, col="black", lwd=lwds)

plot(apclcebu_geo[1:100], apclcebu_fst[1:100], pch=pchs, xlab="Geographic Distance (km)", ylab="", main="APCL Cebu", ylim=ylim, xlim = xlims)
l = lm(apclcebu_fst[lower.tri(apclcebu_fst)] ~ apclcebu_geo[lower.tri(apclcebu_geo)])
lines(apclcebu_geo[1:100][!is.na(apclcebu_geo[1:100])], l$fitted.values, col="black", lwd=lwds)

mod_prbi <- mantel(cebu_geo, cebu_fst, method="pearson", permutations = 10000)
summary(mod<-lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)]))
par(mfrow=c(2,3))
plot(mod, which=1:6)
#mantel(as.dist(cebu_fst) ~ as.dist(cebu_geo), nperm=10000)

mod_apcl <- mantel(apclcebu_geo, apclcebu_fst, method="pearson", permutations = 10000)
mod_apcl
summary(mod<-lm(apclcebu_fst[lower.tri(apclcebu_fst)] ~ apclcebu_geo[lower.tri(apclcebu_geo)]))
par(mfrow=c(2,3))
plot(mod, which=1:6)
#mantel(as.dist(leyte_fst) ~ as.dist(leyte_geo), nperm=10000)



# PRBI vs. APCL fsts: Mantel comparisons (Cebu)

plot(as.dist(apclcebu_fsttrim), as.dist(cebu_fstno15), xlab="APCL Fst/(1-Fst)", ylab="PRBI Fst/(1-Fst)", pch=16, main="APCL vs. PRBI Genetic Distance (no 15)")
mantel(apclcebu_fsttrim, cebu_fstno15)



# PRBI vs. APCL fsts: Mantel comparisons (all data)

plot(as.dist(apclfsttrim), as.dist(fstlintrim), xlab="APCL Fst/(1-Fst)", ylab="PRBI Fst/(1-Fst)", pch=16, main="APCL vs. PRBI Genetic Distance")
mantel(apclfsttrim, fstlintrim)
#abline(coef=c(0,1), col="grey", lty=2)
#l = lm(fst[lower.tri(fst)] ~ apclfst[lower.tri(apclfst)])
#lines(apclfst[lower.tri(apclfst)], l$fitted.values, col="black", lwd=1)


# ANCOVA on PRBI vs. APCL Cebu
x = data.frame(geo=c(apclcebu_geo[lower.tri(apclcebu_fst)], cebu_geo[lower.tri(cebu_fst)]), fst=c(apclcebu_fst[lower.tri(apclcebu_fst)], cebu_fst[lower.tri(cebu_fst)]), spp=c(rep("apcl", length(apclcebu_geo[lower.tri(apclcebu_fst)])), rep("prbi", length(cebu_geo[lower.tri(cebu_fst)]))))
summary(mod<-lm(fst ~ geo*spp, data=x))
plot(x$geo, x$fst, col=c("red", "blue")[as.numeric(x$spp)], pch=16, xlab="Distance (km)", ylab="Fst/(1-Fst)")
abline(prbimod<-lm(x$fst[x$spp=="prbi"]~x$geo[x$spp=="prbi"]),col="blue")
abline(apclmod<-lm(x$fst[x$spp=="apcl"]~x$geo[x$spp=="apcl"]),col="red")
legend("bottomright", legend=unique(x$spp), fill=c("red", "blue"))
summary(apclmod)
summary(prbimod)

summary(mod<-lm(fst ~ geo*spp, data=x))

summary(lm(fst~geo, data=x))





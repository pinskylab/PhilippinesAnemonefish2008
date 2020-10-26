library(vegan)

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/PRBI/090923/")

## Read in data

#fst = read.csv("../090907_noACH_B9/Genepop/PRBI_2009-09-07_fsts_genepop_no12.csv", row.names=1); cebu=8 # Genepop OLD w/out ACH_B9 
#fst = read.csv("../090923/Genepop/PRBI_2009-09-23_fsts_no12.csv", row.names=1); cebu=8 # Genepop
fst = read.csv("../090923/Genepop/PRBI_2009-09-23_fsts_no1215.csv", row.names=1); cebu=7 # Genepop w/out pop15
#fst = read.csv("../090923_noACH_B9/Genepop/PRBI_2009-09-23_noACH_B9_fsts.csv", row.names=1); cebu=8 # Genepop w/out ACH_B9
geo = read.csv("../090907/Arlequin/PRBI_geo_no12_090918.csv", row.names=1) 

## The APCL data for comparison
apclfst = read.csv("../../APCL/090507/genepop/APCL_2009-05-07_fsts_CebuLeyte.csv", header=TRUE, row.names=1)
apclgeo = read.csv("../../APCL/090507/Aclarkii_2009-05-14 geo.csv", header=TRUE, row.names=1)


## Set up data

# linearize fst, remove diagonals
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)

geo[upper.tri(geo, diag=T)] = NA

cebu_fst = as.matrix(fstlin[1:cebu,1:cebu])
cebu_geo = as.matrix(geo[1:cebu,1:cebu])

leyte_fst = as.matrix(fstlin[(cebu+1):(cebu+2),(cebu+1):(cebu+2)])
leyte_geo = as.matrix(geo[9:10, 9:10])

# for APCL
apclfst[upper.tri(apclfst, diag=T)] = NA
apclfstlin = apclfst/(1-apclfst)

apclgeo[upper.tri(apclgeo, diag=T)] = NA
i = is.na(apclgeo)
apclfst[i] = NA

apclcebu_fst = as.matrix(apclfstlin[1:10,1:10])
apclcebu_geo = as.matrix(apclgeo[1:10,1:10])

apclleyte_fst = as.matrix(apclfstlin[11:18,11:18])
apclleyte_geo = as.matrix(apclgeo[11:18, 11:18])


#### Analysis

# Cebu and Leyte mantel together (have to load function below to handle NAs)
par(cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
xlims = c(0,250)
ylim = c(min(leyte_fst, cebu_fst, na.rm=T), max(leyte_fst, cebu_fst, na.rm=T))
col1 = "black" # Cebu
col2 = "grey60" # Leyte
cols = c(rep(col1, 100), rep(col2, 224))

plot(as.dist(geo), as.dist(fst), xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Isolation-by-distance", ylim=ylim, xlim = xlims, pch = 20, lwd=2, col=cols)
x = as.dist(geo)
y = as.dist(fst)
l=(lm(y ~ x))
lines(as.dist(geo)[!is.na(as.dist(geo))], l$fitted.values, col="blue")

summary(l)
predict(l, new = data.frame(x=c(250, 500, 1000)))
mantel.me(fst, geo, nperm=10000, use="pairwise.complete.obs")


# Cebu and Leyte separately
quartz(width=9, height=5)
#par(mfrow=c(1,2)) # for paper?
par(mfrow=c(1,2), cex=1.6, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.9) # for ppt
pchs = 20
lwds = 2
xlims = c(0,250)
ylim = c(min(c(leyte_fst, cebu_fst), na.rm=T), max(c(leyte_fst, cebu_fst), na.rm=T))
plot(cebu_geo[1:64], cebu_fst[1:64], pch=pchs, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Cebu", ylim=ylim, xlim = xlims)
l = lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)])
lines(cebu_geo[1:64][!is.na(cebu_geo[1:64])], l$fitted.values, col="black", lwd=lwds)

plot(leyte_geo[1:64], leyte_fst[1:64], pch=pchs, xlab="Geographic Distance (km)", ylab="", main="Leyte", ylim=ylim, xlim = xlims)
l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)])
lines(leyte_geo[1:64][!is.na(leyte_geo[1:64])], l$fitted.values, col="black", lwd=lwds)

mod_cebu <- mantel(cebu_geo, cebu_fst, method="pearson", permutations = 10000)
mod_cebu
summary(mod<-lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)]))
par(mfrow=c(2,3))
plot(mod, which=1:6)
#mantel(as.dist(cebu_fst) ~ as.dist(cebu_geo), nperm=10000)




# PRBI vs. APCL: Cebu
quartz(width=9, height=5)
#par(mfrow=c(1,2)) # for paper?
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





#########################
## homegrown mantel test to work with missing values
## After Piepho, HP 2005 Permutation tests for the correlation among genetic distances and measures of heterosis. Theor Appl Gen 111:95-99

mantel.me = function(fst, geo, nperm = 1000, use = "pairwise.complete.obs"){
	fstdist = as.dist(fst)
	fstsym = as.matrix(fstdist) # convert to symmetic for permuting
	r = cor(fst[lower.tri(fst)], geo[lower.tri(geo)], use="pairwise.complete.obs") # observed r
	
	n = dim(fst)[1]
	rp = numeric(0)
	for(i in 1:nperm){
		bootn = sample(1:n, n, replace=FALSE)
		bootfst = fstsym[bootn,bootn] # permut rows and columns
		rp = c(rp, cor(bootfst[lower.tri(fst)], geo[lower.tri(geo)], use=use))
	}
	p = (sum(abs(rp)>=abs(r))+1)/(nperm+1)
	out = list(r = r, p=p)
	return(out)
}


####################################
## Fst vs. density
####################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090507/")
fst = read.csv("Aclarkii_2009-05-08 fst.csv", row.names=1)
dens = read.csv("../../../Surveys/surveys2009-03-31.density.csv")
colnames(fst) = row.names(fst)
geo = read.csv("Aclarkii_2009-05-14 geo.csv", row.names=1)
colnames(geo) = row.names(geo)
fst[is.na(geo)] =NA
fstlin = fst/(1-fst)
pops = data.frame(SiteNum=row.names(fst))

# using random site
rand = dens[dens$RandomSite==1,c("SiteNum", "densAPCL")]
rand = merge(rand,pops, all.y=TRUE)
randcol = matrix(rep(rand$densAPCL, length(pops$SiteNum)), byrow=FALSE, nrow=length(pops$SiteNum))
colnames(randcol) = pops$SiteNum
row.names(randcol) = pops$SiteNum
randrow = matrix(rep(rand$densAPCL, length(pops$SiteNum)), byrow=TRUE, nrow=length(pops$SiteNum))
colnames(randrow) = pops$SiteNum
row.names(randrow) = pops$SiteNum
randmean = (randcol+randrow)/2
randmeanlin = randmean[1:(18*18)]
fstlin = as.matrix(fst)[1:(18*18)]
fstlin[fstlin < 0] = 0
i = !is.na(fstlin) & !is.na(randmeanlin)
l = lm(fstlin[i] ~ randmeanlin[i])
plot(randmeanlin[i], fstlin[i], xlab="Random Site Density (Fish/m2)", ylab="Fst")
#plot(randmeanlin[i], fstlin[i]/(1- fstlin[i]), xlab="Random Site Density (Fish/m2)", ylab="Fst/(1-Fst)")
lines(randmeanlin[i], l$fitted.values, col="blue", type="l")
summary(l)

mantel.me(randmean, fst)


# using non-random sites
nr = dens[dens$RandomSite!=1,c("SiteNum", "densAPCL")]
nr = aggregate(nr$densAPCL, by=list(SiteNum=nr$SiteNum), mean, na.rm=TRUE)
nr = merge(nr,pops, all.y=TRUE)
nrcol = matrix(rep(nr$x, length(pops$SiteNum)), byrow=FALSE, nrow=length(pops$SiteNum))
colnames(nrcol) = pops$SiteNum
row.names(nrcol) = pops$SiteNum
nrrow = matrix(rep(nr$x, length(pops$SiteNum)), byrow=TRUE, nrow=length(pops$SiteNum))
colnames(nrrow) = pops$SiteNum
row.names(nrrow) = pops$SiteNum
nrmean = (nrcol+nrrow)/2
nrmeanlin = nrmean[1:(18*18)]
fstlin = as.matrix(fst)[1:(18*18)]
i = !is.na(fstlin)
l = lm(fstlin[i] ~ nrmeanlin[i])
plot(nrmeanlin[i], fstlin[i], , xlab="Non-random Site Density (Fish/m2)", ylab="Fst")
lines(nrmeanlin[i], l$fitted.values, col="blue", type="l")
summary(l)

mantel.me(nrmean, fst)


# Fst vs. Density & Distance: Random Sites
rand = dens[dens$RandomSite==1,c("SiteNum", "densAPCL")]
rand = merge(rand,pops, all.y=TRUE)
randcol = matrix(rep(rand$densAPCL, length(pops$SiteNum)), byrow=FALSE, nrow=length(pops$SiteNum))
colnames(randcol) = pops$SiteNum
row.names(randcol) = pops$SiteNum
randrow = matrix(rep(rand$densAPCL, length(pops$SiteNum)), byrow=TRUE, nrow=length(pops$SiteNum))
colnames(randrow) = pops$SiteNum
row.names(randrow) = pops$SiteNum
randmean = (randcol+randrow)/2
randmeanlin = randmean[1:(18*18)]
fstlin = as.matrix(fst)[1:(18*18)]
geolin = as.matrix(geo)[1:(18*18)]
l = lm(fstlin[i] ~ randmeanlin[i] + geolin[i])
l2 = lm(fstlin[i] ~ geolin[i])
plot(randmeanlin[i], fstlin[i], xlab="Random Site Density (Fish/m2)", ylab="Fst")
lines(randmeanlin[i], l$fitted.values, col="blue", type="l")
summary(l)
summary(l2)

# Fst vs. Density & Distance: Random Sites: Leyte
randmeanlinleyte = randmean[11:18,11:18]
randmeanlinleyte = as.numeric(randmeanlinleyte)
fstlinleyte = as.matrix(fstlin)[11:18,11:18]
fstlinleyte = as.numeric(fstlinleyte)
geolinleyte = as.matrix(geo)[11:18,11:18]
geolinleyte = as.numeric(geolinleyte)
i = !is.na(fstlinleyte) & !is.na(randmeanlinleyte) # not many sites left because don't have site 18
l = lm(fstlinleyte[i] ~ randmeanlinleyte[i] + geolinleyte[i])
i = !is.na(fstlinleyte)
l2 = lm(fstlinleyte[i] ~ geolinleyte[i]) 
summary(l)
summary(l2)
plot(geolinleyte, fstlinleyte)
mantel(geo[11:18, 11:18], fstlin[11:18,11:18])

# Fst vs. Density & Distance: Random Sites: Cebu
randmeanlincebu = randmean[1:10,1:10]
randmeanlincebu = as.numeric(randmeanlincebu)
fstlincebu = as.matrix(fstlin)[1:10,1:10]
fstlincebu = as.numeric(fstlincebu)
geolincebu = as.matrix(geo)[1:10,1:10]
geolincebu = as.numeric(geolincebu)
i = !is.na(fstlincebu) & !is.na(randmeanlincebu)
l = lm(fstlincebu[i] ~ randmeanlincebu[i] + geolincebu[i])
i = !is.na(fstlincebu)
l2 = lm(fstlincebu[i] ~ geolincebu[i])
l2 = lm(fstlincebu ~ geolincebu)
summary(l)
summary(l2)
plot(geolincebu, fstlincebu)
plot(randmeanlincebu, fstlincebu)
plot(randmeanlincebu[i], l2$residuals) # run after trimming to not nas for fstlincebu, randmeanlincebu, fitting l2 with [i]
mantel(geo[1:10, 1:10], fstlin[1:10,1:10])
mantel.me(randmean, fst)



####################################
#### Jackknife over populations
####################################
# See "Arlequin jackknife/" for jackknife over loci

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090507/")

fst = read.csv("Aclarkii_2009-05-08 fst.csv", row.names=1)
#geo = read.csv("Aclarkii_2009-04-24 geo.csv", row.names=1)
geo = read.csv("Aclarkii_2009-05-14 geo.csv", row.names=1)

# linearize fst, remove diagonals
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)
geo[upper.tri(geo, diag=T)] = NA
i = is.na(geo)
fst[i] = NA

cebu_fst = as.matrix(fstlin[1:10,1:10])
cebu_geo = as.matrix(geo[1:10,1:10])
leyte_fst = as.matrix(fstlin[11:18,11:18])
leyte_geo = as.matrix(geo[11:18, 11:18])

library(vegan)

maxpops = 10
mant = data.frame(Leyte_b = numeric(maxpops), Leyte_r = numeric(maxpops), Leyte_r2 = numeric(maxpops), 
Leyte_p = numeric(maxpops), Cebu_b = numeric(maxpops), Cebu_r = numeric(maxpops), Cebu_r2 = numeric(maxpops), Cebu_p = numeric(maxpops))

## Plots and regression
# Leyte
quartz(width=6, height=8)
par(mfrow=c(4,2))
pops = 1:8
popnames = c(18,19,20,22,23,24,25,27)
for(i in pops){
	jackpops = pops[-i]
	jackgen = leyte_fst[jackpops, jackpops]
	jackgeo = leyte_geo[jackpops, jackpops]

	plot(jackgeo[1:49], jackgen[1:49], pch=16, main=paste("W/out", popnames[i]), xlab="Distance (km)", ylab="Fst/(1-Fst)")
	l = lm(jackgen[lower.tri(jackgen)] ~ jackgeo[lower.tri(jackgeo)])
	lines(jackgeo[lower.tri(jackgen)], l$fitted.values, col="orange")

	mant$Leyte_b[i] = l$coefficients[2]
	mant$Leyte_r2[i] = summary(l)$r.squared	
}
mant$Leyte_b[9:10] = NA
mant$Leyte_r[9:10] = NA
mant$Leyte_r2[9:10] = NA
mant$Leyte_p[9:10] = NA

# Cebu
quartz(width=6, height=8)
par(mfrow=c(4,3))
pops = 1:10
popnames = c(7,8,9,10,11,13,14,15,16,17)
for(i in pops){
	jackpops = pops[-i]
	jackgen = cebu_fst[jackpops, jackpops]
	jackgeo = cebu_geo[jackpops, jackpops]

	plot(jackgeo[1:81], jackgen[1:81], pch=16, main=paste("W/out", popnames[i]), xlab="Distance (km)", ylab="Fst/(1-Fst)")
	l = lm(jackgen[lower.tri(jackgen)] ~ jackgeo[lower.tri(jackgeo)])
	lines(jackgeo[lower.tri(jackgen)], l$fitted.values, col="orange")

	mant$Cebu_b[i] = l$coefficients[2]
	mant$Cebu_r2[i] = summary(l)$r.squared	
}



## Mantels
# Leyte
pops = 1:8
for(i in pops){
	jackpops = pops[-i]
	jackgen = leyte_fst[jackpops, jackpops]
	jackgeo = leyte_geo[jackpops, jackpops]

	l = mantel(jackgeo, jackgen, method="pearson", permutations = 10000)
	mant$Leyte_p[i] = l$signif
	mant$Leyte_r[i] = l$statistic
}
# Cebu
pops = 1:10
for(i in pops){
	jackpops = pops[-i]
	jackgen = cebu_fst[jackpops, jackpops]
	jackgeo = cebu_geo[jackpops, jackpops]

	l = mantel(jackgeo, jackgen, method="pearson", permutations = 10000)
	mant$Cebu_p[i] = l$signif
	mant$Cebu_r[i] = l$statistic
}

row.names(mant) = paste("No", 1:10)

# Jackknife analysis on Leyte: See Sokal & Rolf
hist(mant$Leyte_b, breaks=10, col="grey", xlim=c(0,0.0002), freq=F, main="Jackknife estimates of Leyte b", xlab="b")
lines(density(mant$Leyte_b[1:8], adjust=2))

St = summary(mod<-lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)]))$coefficients[2,1]
pops = 1:8
len = length(pops)
ps = numeric(len)
for(i in pops){
	ps[i] = len*St-(len-1)*mant$Leyte_b[i] # the pseudovalues
}
Sthat = mean(ps) # the jackknifed mean
sst = sqrt(sum((St-ps)^2)/(len*(len-1)))

ts = (St-0)/sst # the t-value
dt(ts, df=(len-1)) # the p-value: a mantel test is more appropriate than this, right?

# Jackknife analysis on Cebu: See Sokal & Rolf
hist(mant$Cebu_b, breaks=15, col="grey", xlim=c(0,0.00004), freq=F, main="Jackknife estimates of Cebu b", xlab="b")
lines(density(mant$Cebu_b[1:10], adjust=2)) # adjust = 2x the normal bandwidth

St = summary(mod<-lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)]))$coefficients[2,1]
pops = 1:10
len = length(pops)
ps = numeric(len)
for(i in pops){
	ps[i] = len*St-(len-1)*mant$Cebu_b[i] # the pseudovalues
}
Sthat = mean(ps) # the jackknifed mean
sst = sqrt(sum((St-ps)^2)/(len*(len-1)))

ts = (St-0)/sst # the t-value
dt(ts, df=(len-1)) # the p-value: a mantel test is more appropriate than this, right?


write.csv(mant, paste("MantelPops",Sys.Date(),".csv", sep=""))



#################################################################
## Resampling approach to estimate 95% CI on dispersal distance
#################################################################
library(locfit)

# From Rousset method: slope = 1/4Dsigma2
len = 10000 # how many iterations?

# Island Mean
	m = rnorm(len, mean=3.009e-05, sd = 1.441e-5) # slope of IBD using Cebu & Leyte together (neg Fsts remain), removed APCL285,286 and fixed APCL250, using 090514 GEarth distances

	D = runif(len, 0.16, 252) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0, 252) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0.17, 252) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, He=0.615, stepwise model), to Cebu-Leyte adult density (from Araki et al 2006)
#	D = runif(len, 0.17, 125) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, He=0.615, stepwise model), to 50% Cebu-Leyte adult density (from Araki et al 2006)
#	D = runif(len, 0.06, 252) # spp density from He of Neighborhoods (mu = 1e-3, 150000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0.08, 252) # spp density from He of Neighborhoods (mu = 1e-3, 24000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0.49, 252) # spp density from He of Neighborhoods (mu = 1e-4, 24000km, Theta=2.87, max=100km), to mean Cebu-Leyte adult density
#	D = runif(len, 0.81, 252) # spp density from He of Neighborhoods (mu = 1e-4, 24000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0.16, 400) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=2.87, to mean Cebu-Leyte adult density
#	D = runif(len, 0.0008, 252)

# Leyte
#	m = rnorm(len, mean=9.16e-05, sd=0.000027) # Leyte slope w/ neg Fsts
#	D = runif(len, 0.16, 144) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=3, to mean Leyte adult density

# Cebu
#	m = rnorm(len, mean = 1.910e-05, sd = 1.258e-05) # Ceub slope w/ neg Fsts
#	D = runif(len, 0.16, 317) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, theta=3), to mean Cebu adult density

		
#	D = 10^runif(len, log10(0.16), log10(252)) # spp density from He of Neighborhoods (mu = 5e-4, 24000km, Theta=3, to mean Cebu-Leyte adult density


# Use bounds method for De (upper and lower)
sigma_bounds = sqrt(1/(4*D*m))
	length(sigma_bounds)
	i = !is.na(sigma_bounds) & is.finite(sigma_bounds)
	D = D[i]
	m = m[i]
	sigma_bounds = sigma_bounds[i]
	length(sigma_bounds)

locfit(~sigma_bounds)->fit.sigma
#plot(fit.sigma, xlim=c(0,50))
sigma_bounds[which.max(predict(fit.sigma,newdata=sigma_bounds))]

quantile(sigma_bounds, c(0.025, 0.975))
print(paste("Upper density:", sqrt(1/(4*max(D)*mean(m)))))
print(paste("Lower density:", sqrt(1/(4*min(D)*mean(m)))))
print(paste("Upper density (1000):", sqrt(1/(4*1000*3.009e-05))))
print(paste("Mid density (126):", sqrt(1/(4*126*3.009e-05))))
print(paste("Upper density (252):", sqrt(1/(4*252*3.009e-05))))
print(paste("Lower density (0.16):", sqrt(1/(4*0.16*3.009e-05))))
print(paste("Lower density (0.5):", sqrt(1/(4*0.5*3.009e-05))))


# Use point estimate of De from MNe
mne = 10 # using all pops as source
sigma_point = sqrt(1/(4*mne*m))
	length(sigma_point)
	i = !is.na(sigma_point) & is.finite(sigma_point)
	D = D[i]
	m = m[i]
	sigma_point = sigma_point[i]
	length(sigma_point)

locfit(~sigma_point)->fit.sigma
#plot(fit.sigma, xlim=c(0,50))
sigma_point[which.max(predict(fit.sigma,newdata=sigma_point))]
quantile(sigma_point, c(0.025, 0.975))

# Use SECOND point estimate of De from MNe
mne = 54 # using flanking pops as source
sigma_point2 = sqrt(1/(4*mne*m))
	length(sigma_point)
	i = !is.na(sigma_point) & is.finite(sigma_point)
	D = D[i]
	m = m[i]
	sigma_point = sigma_point[i]
	length(sigma_point)

locfit(~sigma_point)->fit.sigma
#plot(fit.sigma, xlim=c(0,50))
sigma_point[which.max(predict(fit.sigma,newdata=sigma_point))]
quantile(sigma_point, c(0.025, 0.975))

print(paste("Density = 10 (Mne):", sqrt(1/(4*10*3.009e-05))))
print(paste("Density = 54 (Mne):", sqrt(1/(4*54*3.009e-05))))
print(paste("Density = 1000 (absolute):", sqrt(1/(4*1000*3.009e-05))))
print(paste("Kinlan & Gaines 2003:", 0.0016*3.009e-5^-1.0001))

#par(mfrow=c(4,1))
par(mfrow=c(3,1), cex=1.2, cex.axis=0.7, omi = c(0,0,0,0), mgp=c(1.9,0.7,0), mar=c(4,3,1.5,2))
plot(-10,-10, xlim=c(-5,1.1*max(D)), ylim=c(0, 1.5/(max(D)-min(D))), ylab="Density", main="Effective density", xlab="Adults per km", bty="L", yaxp=c(0,0.006, 2))
rect(min(D), 0, max(D), 1/(max(D)-min(D)), col="grey40", border="NA")
x = seq(min(m),max(m),by=(max(m)-min(m))/100)
y = 1/(sd(m)*sqrt(2*pi))*exp(-(x-mean(m))^2/(2*sd(m)^2)) # plot a normal distribution
plot(x,y,type="l", col=NA, bty="L", main="Slope", xlab="b", ylab="Density", yaxp=c(0,30000,2))
polygon(x,y, col="grey40", border=NA)
hist(sigma, breaks=c(seq(0,200, by=2), max(sigma)), xlab = "Spread (km) (sigma)", xlim=c(0,100), freq=FALSE, col="grey40", main="Spread", yaxp=c(0,0.12, 2))

# Plot bounds and point estimate methods together as histograms
maxspr = max(c(sigma_bounds, sigma_point))
xlim = c(0,200)
xlim = c(0,50)
log = "x"
by =2
by =1
main = "Bounding method for De"
main = "Cebu w/ mean density"
main = "Leyte w/ mean density"
main = "Cebu w/ Cebu density"
main = "Leyte w/ Leyte density"
#quartz(width=12, height=5)
#par(mfrow=c(2,1))
quartz(width=12, height=8) # for plotting 4
par(mfrow=c(4,1))
par(cex=1.8, cex.axis=0.7, omi = c(0,0,0,0), mgp=c(1.4,0.4,0), mar=c(2.5,3,1.5,2))
hist(sigma_bounds, breaks=c(seq(0,200, by=by), maxspr), xlab = "", ylab="", xlim=xlim, freq=FALSE, col="grey20", border="grey20", main= main, yaxp=c(0,0.12, 2))
par(mar=c(2.5,3,1.5,2))
hist(sigma_point, breaks=c(seq(0,200, by=2), maxspr), xlab = "Spread (km) (sigma)", ylab="", xlim=xlim, freq=FALSE, col="grey20", border="grey20", main="Point estimate for De", yaxp=c(0,0.12, 2))

# Plot bounds and point method on log scale
maxspr = max(c(sigma_bounds, sigma_point))
ticksm = c(seq(0.1,0.9,by=0.1), seq(1,9, by=1), seq(10,90, by=10), seq(100, 1000, by=100))
ticklg = c(0.1,1,10,100,1000)
quartz(width=12, height=6)
par(mfrow=c(3,1))
par(cex=1.8, cex.axis=0.8, omi = c(0,0,0,0), mgp=c(1.4,0.4,0), mar=c(2.5,3,0,2))
hist(log10(sigma_bounds), breaks=c(seq(-1,log10(200), by=0.05), log10(maxspr)), xlab="", main="", axes=FALSE, ylab="", col="grey20", border="grey20")
axis(1, labels=F, at=log10(ticksm), tcl=-0.3)
axis(1, labels=ticklg, at=log10(ticklg))
hist(log10(sigma_point), breaks=c(seq(-1,log10(200), by=0.05), log10(maxspr)), xlab="", main="", axes=FALSE, ylab="", col="grey20", border="grey20")
axis(1, labels=F, at=log10(ticksm), tcl=-0.3)
axis(1, labels=ticklg, at=log10(ticklg))
hist(log10(sigma_point2), breaks=c(seq(-1,log10(200), by=0.05), log10(maxspr)), xlab="", main="", axes=FALSE, ylab="", col="grey20", border="grey20")
axis(1, labels=F, at=log10(ticksm), tcl=-0.3)
axis(1, labels=ticklg, at=log10(ticklg))


require(lattice)
contourplot(sigma ~ D+m, region=TRUE, contour=TRUE)




########### Compare mean and standard deviation
# Lockwood defined mean as mean of half the distribution (e.g., the positive half)

# sd = 10

# normal
u = rnorm(10000, 0, 7)
sd = sd(u)
mean(u[u>0])
mean(u[u>0])/sd # 80%

# lognormal
u = rlnorm(10000, 0, 1.526)
v = -rlnorm(10000, 0, 1.526)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 33%

# laplace
sd = 10
u = runif(len, -1/2, 1/2)
u = -sd/sqrt(2)*sign(u)*log(1-2*abs(u))
#hist(u)
sd(u)
mean(u[u>0])
mean(u[u>0])/sd(u) # 70%

# Weibull
u = rweibull(10000, .408)
v = -rweibull(10000, .408)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 32%

# Gamma
u = rgamma(10000, 9.6)
v = -rgamma(10000, 9.6)
sd(c(u,v))
mean(u[u>0])
mean(u[u>0])/sd(c(u,v)) # 95%

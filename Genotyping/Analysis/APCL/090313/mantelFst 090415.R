#setwd("C:/Documents and Settings/Mollie Manier/Desktop/Users/Malin/Philippines 2008/Analysis/090313/")
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")

fst = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 fst.csv", row.names=1)
geo = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 geo.csv", row.names=1)

# linearize fst, remove diagonals, (maybe: and set negs to zeros)
#fst[fst < 0 & !is.na(fst)] = 0
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
#library(ecodist)
#library(ape)

# Cebu and Leyte mantel together (have to load function below to handle NAs)
mantel.me(fst, geo, nperm=10000, use="pairwise.complete.obs")
plot(as.dist(geo), as.dist(fst))
l=(lm(as.dist(fst) ~ as.dist(geo)))
lines(as.dist(geo)[!is.na(as.dist(geo))], l$fitted.values, col="blue")
summary(l)


# Cebu and Leyte separately
quartz(width=9, height=5)
par(mfrow=c(1,2))
ylim = c(min(leyte_fst, cebu_fst, na.rm=T), max(leyte_fst, cebu_fst, na.rm=T))
plot(cebu_geo[1:100], cebu_fst[1:100], pch=1, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Cebu", ylim=ylim)
l = lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)])
lines(cebu_geo[1:100][!is.na(cebu_geo[1:100])], l$fitted.values, col="blue")

mantel(cebu_geo, cebu_fst, method="pearson", permutations = 10000)
summary(lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)]))$coefficients
#mantel(as.dist(cebu_fst) ~ as.dist(cebu_geo), nperm=10000)

plot(leyte_geo[1:64], leyte_fst[1:64], pch=1, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Leyte", ylim=ylim)
l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)])
lines(leyte_geo[1:64][!is.na(leyte_geo[1:64])], l$fitted.values, col="blue")

mantel(leyte_geo, leyte_fst, method="pearson", permutations = 10000)
summary(lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)]))$coefficients
#mantel(as.dist(leyte_fst) ~ as.dist(leyte_geo), nperm=10000)



# remove pop #18 from Leyte: plot Leyte with Cebu
quartz(width=9, height=5)
par(mfrow=c(1,2))
ylims =c(-0.015, 0.03)
plot(cebu_geo[1:100], cebu_fst[1:100], pch=1, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Cebu", ylim=ylims)
l = lm(cebu_fst[lower.tri(cebu_fst)] ~ cebu_geo[lower.tri(cebu_geo)])
lines(cebu_geo[1:100][!is.na(cebu_geo[1:100])], l$fitted.values, col="blue")
summary(l)

leyte_fstno18 = as.matrix(fstlin[12:18,12:18])
leyte_geono18 = as.matrix(geo[12:18, 12:18])

plot(leyte_geo[1:64], leyte_fst[1:64], pch=1, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="Leyte", ylim=ylims)
l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)])
lines(leyte_geo[1:64][!is.na(leyte_geo[1:64])], l$fitted.values, col="blue")
summary(l)

points(leyte_geono18[1:49], leyte_fstno18[1:49], pch=16)
l = lm(leyte_fstno18[lower.tri(leyte_fstno18)] ~ leyte_geono18[lower.tri(leyte_geono18)])
lines(leyte_geono18[1:49][!is.na(leyte_geono18[1:49])], l$fitted.values, col="orange")
summary(l)

mantel(cebu_geo, cebu_fst, method="pearson", permutations=10000)
mantel(leyte_geo, leyte_fst, method="pearson", permutations=10000)
mantel(leyte_geono18, leyte_fstno18, method="pearson", permutations = 10000)


# use a longer distance measure for Pop 18 (Pintuyan): 100km instead of 25km
leyte_geolong = leyte_geo + cbind(rep(75,8), matrix(0,nrow=8,ncol=7))

mantel(leyte_geolong, leyte_fst, method="pearson", permutations = 10000)
summary(lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geolong[lower.tri(leyte_geolong)]))

plot(leyte_geo[1:64], leyte_fst[1:64], pch=1, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", xlim=c(0,300), main="Leyte")
l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geo[lower.tri(leyte_geo)])
lines(leyte_geo[lower.tri(leyte_geo)], l$fitted.values, col="blue")
summary(l)
points(leyte_geolong[lower.tri(leyte_geolong)],leyte_fst[lower.tri(leyte_fst)], pch=16)
l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geolong[lower.tri(leyte_geolong)])
lines(leyte_geolong[lower.tri(leyte_geolong)], l$fitted.values, col="orange")
summary(l)



################################################
# convert b to dispersal distance
# using Leyte w/ 100km Pintuyan-Padre Burgos

l = lm(leyte_fst[lower.tri(leyte_fst)] ~ leyte_geolong[lower.tri(leyte_geolong)])
i = summary(l)$coefficients

dens1 = 143/46.6/600 # fish/km, from lowest confidence interval of MLNE estimate for Philippines (Ne=143, 3/17), then divide by highest Philippines to CebuLeyte ratio from GIS area calculations (46.6), then by 600 km for Cebu/Leyte coastlines
#dens2 = 135.5 # fish/km, from processPhils2009-03-26.R, using polygon areas and matching APCL density surveys to arcs or polygons, times 10% Ne/N ratio (Frankham 1995 Gen Res)
dens2 = 235 # fish/km from processPhils2009-03-26.R, converting mean adults/m2 across Random Sites to adults/km with 150m reef width
dens3 = 1000/42.2/600 # approximate average across MLNE and Waples 1989 methods, and ave Philippines/CebuLeyte ratio (42.2)
#dens4 = 13.5 # like dens2, but using 10e-2 Ne/N ratio (more appropriate for marine spp? 10e-2 to 10e-6?)
dens4 = 235/1000 # like dens2, but using an Ne/N ratio

dispersal = data.frame(name = c("mean", "UL", "LL"), b = c(i[2,1], i[2,1]+1.96*i[2,2], i[2,1]-1.96*i[2,2]))
dispersal$binv = 1/dispersal$b
dispersal$sigma_Ne = sqrt(dispersal$binv/4/dens1) # using the smallest effective density estimate
dispersal$sigma_Ecol = sqrt(dispersal$binv/4/dens2) # using the biggest ecological density estimate with Ne/N = 10%
dispersal$sigma_Ne2 = sqrt(dispersal$binv/4/dens3) # using the approx genetic Ne estimate (1000)
dispersal$sigma_Ecol2 = sqrt(dispersal$binv/4/dens4) # using a different ecological density estimate
format(dispersal, digits=2)



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
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")
fst = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 fst.csv", row.names=1)
dens = read.csv("../../../Surveys/surveys2009-03-31.density.csv")
colnames(fst) = row.names(fst)
geo = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 geo.csv", row.names=1)
colnames(geo) = row.names(geo)
fst[is.na(geo)] =NA
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
i = !is.na(fstlin) & !is.na(randmeanlin)
l = lm(fstlin[i] ~ randmeanlin[i])
plot(randmeanlin[i], fstlin[i], xlab="Random Site Density (Fish/m2)", ylab="Fst")
#plot(randmeanlin[i], fstlin[i]/(1- fstlin[i]), xlab="Random Site Density (Fish/m2)", ylab="Fst/(1-Fst)")
lines(randmeanlin[i], l$fitted.values, col="blue", type="l")
summary(l)

mantel.me(randmean, fst)




####################################
#### Jackknife over populations
####################################
# See "Arlequin jackknife/" for jackknife over loci

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")

fst = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 fst.csv", row.names=1)
geo = read.csv("Aclarkii_2009-03-13 CebuLeyte noACH_A7 geo.csv", row.names=1)

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

write.csv(mant, paste("MantelPops",Sys.Date(),".csv", sep=""))




#################################################################
## Resampling approach to estimate 95% CI on dispersal distance
#################################################################

# From Rousset method: slope = 1/4Dsigma2
len = 10000 # how many iterations?

#	m = rnorm(len, mean=6.90e-05, sd=0.000024) # slope of IBD using Leyte w/ pop 18 (neg Fsts to zero)
#	m = rnorm(len, mean=2.56e-05, sd=0.00001107) # slope of IBD using Cebu and Leyte together (neg Fsts to zero)
	m = rnorm(len, mean=9.16e-05, sd=0.000027) # slope of IBD using Leyte w/ pop 18 (neg Fsts remain)

#	D = runif(len, 5, 144) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Philippines to Leyte adult density
	D = runif(len, 5, 72) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Philippines to 50% of Leyte adult density
#	D = runif(len, 1.8, 72) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Philippines to 50% of Leyte adult density
#	D = runif(len, 1.8, 144) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Whole Range to Leyte adult density
#	D = runif(len, 0.01, 14) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6, range=Whole Range, and 2Dimensional, up to 10% of Leyte adult density
#	D = 10^(runif(len,0.699,2.15)) # log10-uniform from 5 to 144
#	D = 10^(runif(len,0.2553, 2.15)) # log10-uniform from 1.8 to 144
#	D = 10^(runif(len, -2 ,2.15)) # log10-uniform from 0.01 to 144
		
	sigma = sqrt(1/(4*D*m))

length(sigma)
i = !is.na(sigma) & is.finite(sigma)
D = D[i]
m = m[i]
sigma = sigma[i]
length(sigma)

as.list(density(sigma))$x[which.max(as.list(density(sigma))$y)]
quantile(sigma, c(0.025, 0.975))
#par(mfrow=c(4,1))
par(mfrow=c(3,1))
hist(D, xlab = "Effective density prior (Adults/km)", main="Density")
#plot(density(log10(D)), xlab = "Effective density prior (Adults/km)", main="Density")
hist(m, xlab = "Slope", main = "Isolation by Distance Relationship")
#plot(density(sigma), xlab = "Dispersal distance (km) (sigma)", "Dispersal distance")
plot(density(sigma), xlab = "Dispersal distance (km) (sigma)", "Dispersal distance", xlim=c(0,100))

plot(density(sigma), xlim=c(0,100))

require(lattice)
contourplot(sigma ~ D+m, region=TRUE, contour=TRUE)



# From Wright 1969 method: N = kDsigma
len = 1000 # how many iterations?

N = runif(len, 50 , 1000)

k = runif(len, 1.46, 3.545)

#	D = runif(len, 5, 144) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Philippines to Leyte adult density
#	D = runif(len, 1.8, 144) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6 and range=Whole Range to Leyte adult density
#	D = runif(len, 0.01, 14) # spp density, from He of Neighborhoods 090407.R lower limit for He =0.6, range=Whole Range, and 2Dimensional, up to 10% of Leyte adult density
	D = 10^(runif(len,0.699,2.15)) # log10-uniform from 5 to 144
#	D = 10^(runif(len, -2 ,2.15)) # log10-uniform from 0.01 to 144

sigma =N/(k*D)

length(sigma)
sigma = sigma[!is.na(sigma) & is.finite(sigma)]
length(sigma)

as.list(density(sigma))$x[which.max(as.list(density(sigma))$y)]
par(mfrow=c(2,2))
plot(density(sigma), log="x")
plot(density(N))
plot(density(k))
plot(density(D))
plot(density(sigma), xlim=c(0,100))

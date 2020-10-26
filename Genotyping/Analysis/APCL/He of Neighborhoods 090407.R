setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")


# Plot He against Density for various Neighborhood Nes

quartz(width=9, height=5)
par(mfrow=c(1,2))

i = (he$mu == 1e-5 & he$Phils.km==15000)

xlim = c(.01,400)
ylim = c(0,1)
j = he$N.1.[i] == 100
plot(he$N1.km[i][j], he$He.s.[i][j], log="x", xlim=xlim, ylim=ylim, col="black", type="b", xlab="Effective number of adults per km", ylab = "Species Heterozygosity (He)", main="Philippines (15,000 km)")
j = he$N.1.[i] == 1000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="blue")
j = he$N.1.[i] == 10000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="green")
j = he$N.1.[i] == 100000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="orange")

abline(h=0.6, lty = 2)

legend("topleft", legend=c("100", "1,000", "10,000", "100,000"), pch=1, lty=1, col=c("black", "blue", "green", "orange"), title="Neighborhood Ne")



i = (he$mu == 1e-5 & he$Phils.km==150000)

xlim = c(.01,400)
ylim = c(0,1)
j = he$N.1.[i] == 100
plot(he$N1.km[i][j], he$He.s.[i][j], log="x", xlim=xlim, ylim=ylim, col="black", type="b", xlab="Effective number of adults per km", ylab = "Species Heterozygosity (He)", main="Whole Range (150,000 km)")
j = he$N.1.[i] == 1000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="blue")
j = he$N.1.[i] == 10000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="green")
j = he$N.1.[i] == 100000
points(he$N1.km[i][j], he$He.s.[i][j], type="b", col="orange")

abline(h=0.6, lty = 2)

legend("topleft", legend=c("100", "1,000", "10,000", "100,000"), pch=1, lty=1, col=c("black", "blue", "green", "orange"), title="Neighborhood Ne")



########################################
# Plot He vs. Density for various sigma
cols = c("gray0", "gray10", "gray20", "gray30", "gray40", "gray50", "gray60", "darkorange", "darkorange4")

quartz(width=9, height=5)
par(mfrow=c(1,2))

i = which(he$mu == 1e-5 & he$Phils.km==15000)
dists = unique(he$Sigma[i])
dists = dists[!is.na(dists)]
len = length(dists)

xlim = c(.01,400)
ylim = c(0,1)
j = he$Sigma[i] == dists[1]
plot(he$N1.km[i][j], he$He.s.[i][j], log="x", xlim=xlim, ylim=ylim, col=cols[1], type="b", xlab="Effective number of adults per km", ylab = "Species Heterozygosity (He)", main="Philippines (15,000 km)")
for(k in 2:len){
	j = he$Sigma[i] == dists[k]
	print(sum(j, na.rm=T))
	points(he$N1.km[i][j], he$He.s.[i][j], type="b", col=cols[k])
}

abline(h=0.6, lty = 2)
legend("topleft", legend=signif(dists,2), pch=1, lty=1, col=cols[1:len], title="Sigma (km)")


i = which(he$mu == 1e-5 & he$Phils.km==150000)
dists = unique(he$Sigma[i])
dists = dists[!is.na(dists)]
len = length(dists)

xlim = c(.01,400)
ylim = c(0,1)
j = he$Sigma[i] == dists[1]
plot(he$N1.km[i][j], he$He.s.[i][j], log="x", xlim=xlim, ylim=ylim, col=cols[1], type="b", xlab="Effective number of adults per km", ylab = "Species Heterozygosity (He)", main="Whole Range (150,000 km)")
for(k in 2:len){
	j = he$Sigma[i] == dists[k]
	print(sum(j, na.rm=T))
	points(he$N1.km[i][j], he$He.s.[i][j], type="b", col=cols[k])
}

abline(h=0.6, lty = 2)





#########################################################################
### Contour plots of Species-wide He that are prettier: use linear model
library(lattice)

dens = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 200)
sigma = c(.1, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)

km = 150000 # whole spp range
mu = 1e-5

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma*3.545 # Neighborhood length
he$S = km/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns = sqrt(he$S)*he$Nn # spp Ne
he$thetaS = 4*he$Ns*mu
he$He.s = he$thetaS/(1+he$thetaS)

contourplot(He.s ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance (km)", main = "Species Heterozygosity \n (whole range 150,000km, mu=1e-5)", panel = function(...){ 
	panel.contourplot(...)
	panel.abline(3,0, lty=2)
	panel.abline(v = .25, lty=2)
})


# For the Philippines (15k km)
#km = 15000
km = 24000 # kms of reef, from ReefBase GIS
#km = 36000 # CIA World Factbook
#km = 150000 # whole range in SE Asia, from Reefbase GIS
#mu = 1e-5
#mu = 1e-4
mu = 5e-4

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma*3.545 # Neighborhood length (Gaussian kernel)
#he$L = he$sigma*1.461 # Neighborhood length (n = 1/2 in Bateman leptokurtic kernel)
he$S = km/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns = sqrt(he$S)*he$Nn # spp Ne
he$thetaS = 4*he$Ns*mu
he$He.s = he$thetaS/(1+he$thetaS)

contourplot(He.s ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance sigma (km)", main = paste("Species Heterozygosity \n(Philippines ",  km, "km, mu=", mu, ")",sep=""), panel = function(...){ 
	panel.contourplot(...)
	panel.abline(3,0, lty=2)
	#panel.abline(v = .71, lty=2) # for km=15000 and mu = 1e-5
	#panel.abline(v = -0.28, lty=2) # for km=15000 and mu = 1e-4
	#panel.abline(v = -0.98, lty=2) # for km=15000 and mu = 5e-4
	panel.abline(v = -1.08, lty=2) # for km=36000 and mu = 5e-4
	#panel.abline(v = -1.18, lty=2) # for km=36000 and mu = 5e-4
	#panel.abline(v = -1.5, lty=2) # for km=36000 and mu = 5e-4
})

10^(0.71)
10^(-0.28)
10^(-0.98)
10^(-1.08)
10^(-1.18)
10^(-1.5)




###################################################
### Contour plots of Species-wide He: use 2D model
dens = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 200)
sigma = c(.1, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)

km = 900000 # km2 for Phils
mu = 1e-5

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma^2*pi*4 # Neighborhood area p.303 in Wright 1969
he$S = km/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns =he$S*he$Nn # spp Ne, p299 in Wright 1969
he$thetaS = 4*he$Ns*mu
he$He.s = he$thetaS/(1+he$thetaS)

contourplot(He.s ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance (km)", main = "Species Heterozygosity \n (Philippines 900k km2, mu=1e-5)", panel = function(...){ 
	panel.contourplot(...)
	#panel.abline(3,0, lty=2)
	#panel.abline(v = .25, lty=2)
})


# For the whole range 37000000 km2
km = 37000000
mu = 1e-5

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma^2*pi*4 # Neighborhood area p.303 in Wright 1969
he$S = km/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns =he$S*he$Nn # spp Ne, p299 in Wright 1969
he$thetaS = 4*he$Ns*mu
he$He.s = he$thetaS/(1+he$thetaS)

contourplot(He.s ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance sigma (km)", main = "Species Heterozygosity \n(whole range 37M km2, mu=1e-5)", panel = function(...){ 
	panel.contourplot(...)
	#panel.abline(3,0, lty=2)
	#panel.abline(v = .7, lty=2)
	#panel.abline(log10(300),0, lty=2)
	#panel.abline(v = -0.0, lty=2)
})





#########################################################################
### Contour plots of Neighborhood size (Nn)
library(lattice)

dens = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 200)
sigma = c(.1, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
# can probably use expand.grid() in the above line
he$L = he$sigma*3.545 # Neighborhood length
he$Nn = he$dens*he$L # Neighborhood Ne

contourplot(Nn ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,1,10,160,1000,10000), col.regions = heat.colors(5), scales=list(log = 10), xlab="Effective adults per km (De)", ylab = "Dispersal spread (sigma) (km)", main = "Neighborhood Size (Nn)", panel = function(...){ 
	panel.contourplot(...)
	panel.abline(3,0, lty=2)
	panel.abline(v = -1.4, lty=2)
})

10^-1.4


# a different version, not using lattice (doesn't work well)
dens2 = dens
for(i in 2:length(sigma)){
	dens2 = rbind(dens2, dens)	
}
sigma2 = sigma
for(i in 2:length(dens)){
	sigma2 = cbind(sigma2, sigma)	
}
dim(sigma2)
dim(dens2)
Nn = matrix(nrow=length(sigma), ncol = length(dens))
Nn = dens2 * sigma * 3.545
Nn = t(Nn)

filled.contour(log(dens), log(sigma), Nn)
contour(log(dens), log(sigma), Nn, levels=160, add=TRUE)
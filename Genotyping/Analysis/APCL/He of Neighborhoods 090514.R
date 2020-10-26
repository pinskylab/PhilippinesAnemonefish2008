setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")




#########################################################################
### Contour plots of Species-wide He that are prettier: use linear model
## Now use stepwise-mutation model
library(lattice)

dens = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 200)
sigma = c(.1, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)


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
#he$thetaS = 4*he$Ns*mu 
#he$He.s = he$thetaS/(1+he$thetaS) # infinite-alleles model
he$He.s = 1-1/sqrt(8*he$Ns*mu+1) # stepwise model

contourplot(He.s ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1, 0.2, 0.4, 0.615, 0.8, 0.9, 1), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance sigma (km)", main = paste("Species Heterozygosity \n(Philippines ",  km, "km, mu=", mu, ", stepwise model)",sep=""), panel = function(...){ 
	panel.contourplot(...)
	panel.abline(3,0, lty=2)
	panel.abline(v = -0.78, lty=2) # for km=24000 and mu = 5e-4
})

10^(-0.78)





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



#########################################################################
### Contour plots of Species-wide THETA that are prettier: use linear (1D) model
library(lattice)

dens = c(0.001, 0.005, 0.01, 0.05, 0.075,  0.1, 0.2, 0.3, 0.5, 1, 5, 10, 50, 100, 200)
sigma = c(.1, 1, 5, 10, 50, 100, 500, 750, 1000, 2000, 3000, 5000, 10000)

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)


# For the Philippines (15k km)
#km = 15000
km = 24000 # kms of reef, from ReefBase GIS
#km = 36000 # CIA World Factbook
#km = 150000 # whole range in SE Asia, from Reefbase GIS
#mu = 1e-5
#mu = 1e-4
mu = 1e-4
mu = 5e-4
mu = 10e-4

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma*3.545 # Neighborhood length (Gaussian kernel)
#he$L = he$sigma*1.461 # Neighborhood length (n = 1/2 in Bateman leptokurtic kernel)
he$S = km/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns = sqrt(he$S)*he$Nn # spp Ne
he$thetaS = 4*he$Ns*mu 

contourplot(thetaS ~ dens+sigma, data=he, region=TRUE, contour=TRUE, at = c(0,0.1,2.53,3.27,100), col.regions = heat.colors(8), scales=list(log = 10), xlab="Effective adults per km", ylab = "Dispersal distance sigma (km)", main = paste("Species Theta \n(Philippines ",  km, "km, mu=", mu, ")",sep=""), panel = function(...){ 
	panel.contourplot(...)
	panel.abline(3,0, lty=2)
	panel.abline(v = -0.1, lty=2) # for km=24000 and mu = 5e-4
	panel.abline(v = -0.78, lty=2) # for km=24000 and mu = 5e-4
	panel.abline(v = -1.1, lty=2) # for km=24000 and mu = 10e-4
})

10^(-0.1)
10^(-0.78)
10^(-1.1)


#######################################################
## Plot De as function of sigma for a given theta and various mus

sigma = c(.1, 1, 5, 10, 50, 100, 500, 750, 1000, 2000, 3000, 5000, 10000)

k = 24000 # kms of reef in Phils, from Reefbase GIS
mu = c(1e-4, 5e-4, 1e-3) # possible mutation rates
#mu = c(1e-6, 1e-4, 5e-4, 1e-3, 1e-2) # possible mutation rates
#theta = 3 # observed theta
theta = 3.55 # observed theta 5/31/09 (missing 24,25,27)
a = 3.545 # Gaussian dispersal

sigmamax = 1000 # max sigma to consider
ymin = 0.001

if(!(sigmamax %in% sigma)){
	print("You didn't include sigmamax in sigma!!")
}

quartz(width=6, height=6)
col = heat.colors(length(sigma))
plot(.00001,.00001, ylim=c(ymin, 300), xlim=c(0.9*min(sigma), 1.1*max(sigma)), log="xy", ylab = "Effective density, De (individuals/km)", xlab = "Dispersal spread, sigma (km)", main=paste("Effective density in continous populations\n(theta = ", theta, ", k = ", k, ", a = ", a, ")", sep=""))
for(i in 1:length(mu)){
	De = theta/(4*mu[i]*sqrt(a*sigma*k))
	lines(sigma, De, col=col[i])
	Decrit = De[sigma==sigmamax]
	lines(x = c(sigmamax, sigmamax), y = c(ymin, Decrit), lty = 2, col="black")
	lines(x=c(0.9*min(sigma), 1000), y = c(Decrit, Decrit), lty=2)
	print(paste("mu = ", mu[i], ", Decrit = ", Decrit, sep=""))
}

legend("topright", legend = mu, fill=col, title="Mutation rate")

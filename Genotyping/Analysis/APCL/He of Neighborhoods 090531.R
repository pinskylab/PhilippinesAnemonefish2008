setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090313/")




#########################################################################
### Contour plots of Species-wide He that are prettier: use linear model
## Now use stepwise-mutation model
library(lattice)


#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)

#dens = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 200)
#sigma = c(.1, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000)

densexp = seq(-3, 2.3, by=0.05)
dens = 10^densexp
sigmaexp = seq(-1,4, by=0.1)
sigma=10^sigmaexp


# For the Philippines (15k km)
#k = 15000
k = 24000 # kms of reef, from ReefBase GIS
k = 500
#k = 36000 # CIA World Factbook
#k = 150000 # whole range in SE Asia, from Reefbase GIS
#mu = 1e-5
#mu = 1e-4
mu = 5e-4
sigmamax = 1000 # maximum sigma to consider reasonable
He = 0.615 # observed He
a = 3.545 # dispersal kernel shape parameter (3.545 for Gaussian, 1.461 for Bateman leptokurtic with n=1/2)

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma*a # Neighborhood length
he$S = k/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns = sqrt(he$S)*he$Nn # spp Ne
#he$thetaS = 4*he$Ns*mu 
#he$He.s = he$thetaS/(1+he$thetaS) # infinite-alleles model
he$He.s = 1-1/sqrt(8*he$Ns*mu+1) # stepwise model

Decrit = 1/((1-He)^2*8*mu*sqrt(a*sigmamax*k))-1/(8*mu*sqrt(a*sigmamax*k))
Decritexp = log10(Decrit)

#colorbreaks = c(0,0.1, 0.2, 0.4, 0.615, 0.8, 0.9, 1) # few color breaks
colorbreaks = seq(0,1,by=0.005) # nearly continuous
#cols = heat.colors(length(colorbreaks)) # for ppt
cols = gray(seq(from=0.1, to=1, length.out=length(colorbreaks))) # for journal article
lwdlines = 3 # line weight on the graph: ppt
#axislwd=1 # line weight for axis and tickmarks: paper
axislwd =2# ppt
#fontsize = 12 # base font size: paper
fontsize = 22 # ppt
axiscex = 1 # expansion size for axis font: paper
axiscex = 0.75 # expansion size for axis font: ppt
labelcex = 1 # expansion size for label font


#title = paste("Species Heterozygosity \n(Philippines ",  km, "k, mu=", mu, ", stepwise model)",sep="") # main title for notebook
title = "Species Heterozygosity"

trellis.par.set(par.xlab.text=list(cex=labelcex),par.ylab.text=list(cex=labelcex),par.zlab.text=list(cex=labelcex), fontsize=list(text=fontsize))

#quartz(height=7, width=7.5)
contourplot(He.s ~ sigma+dens, data=he, region=TRUE, contour=FALSE, at=colorbreaks, col.regions = rev(cols), scales=list(log = 10, cex=axiscex), ylab="Effective density (indivs/km)", xlab = "Dispersal distance (km)", main = title, par.settings = list(axis.line = list(lwd = axislwd)), aspect="iso", panel = function(at,contour,region,lwd,...){ 
	panel.contourplot(contour=FALSE,at=colorbreaks,...) # plot the underlying colors
	panel.contourplot(at=c(0,He),region=FALSE, contour=TRUE,lwd=lwdlines,...) # add a contour`
	panel.segments(3,min(densexp),3,Decritexp, lty=2, lwd=lwdlines)
	panel.segments(min(sigmaexp), Decritexp,3, Decritexp, lty=2, lwd=lwdlines)
})

# commands to find which lattice settings to modify
show.settings()
str(trellis.par.get())


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

#dens = c(.1, 1, 2, 10, 20, 33.33333333, 100, 200)
#sigma = c(0.28, 0.85, 1.41, 2.82, 14.1, 28.21, 282.09)

#dens = c(0.001, 0.005, 0.01, 0.05, 0.075,  0.1, 0.2, 0.3, 0.5, 1, 5, 10, 50, 100, 200)
#sigma = c(.1, 1, 5, 10, 50, 100, 500, 750, 1000, 2000, 3000, 5000, 10000)

densexp = seq(-3, 2.3, by=0.05)
dens = 10^densexp
sigmaexp = seq(-1,4, by=0.1)
sigma=10^sigmaexp


# For the Philippines (15k km)
#k = 15000
k = 24000 # kms of reef, from ReefBase GIS
#k = 36000 # CIA World Factbook
#k = 150000 # whole range in SE Asia, from Reefbase GIS
#mu = 1e-5
#mu = 1e-4
mu = 5e-4
#mu = 10e-4

#sigmamax = 1000 # maximum sigma to consider reasonable
sigmamax = 100 # maximum sigma to consider reasonable

He = 0.615 # observed heterozygosity
theta = 1/(2*(1-He)^2)-1/2
theta = 2.54

#theta = 3.5 # observed theta
thetalog = log10(theta) # observed theta
a = 3.545 # dispersal kernel shape parameter (3.545 for Gaussian, 1.461 for Bateman leptokurtic with n=1/2)

he = data.frame(dens=rep(dens, length(sigma)), sigma = rep(sigma,rep(length(dens), length(sigma))))
he$L = he$sigma*a # Neighborhood length (Gaussian kernel)
he$S = k/he$L  # Number of neighborhoods
he$Nn = he$dens*he$L # Neighborhood Ne
he$Ns = sqrt(he$S)*he$Nn # spp Ne
he$thetaS = 4*he$Ns*mu 
he$logthetaS = log10(4*he$Ns*mu)

Decrit = theta/(4*mu*sqrt(a*sigmamax*k))
Decritexp = log10(Decrit)

#colorbreaks = c(0,0.1, 0.2, 0.4, 0.615, 0.8, 0.9, 1)
colorbreaksexp = seq(min(log10(he$thetaS)),max(log10(he$thetaS)),by=0.05)
#colorbreaks = 10^colorbreaksexp
colorbreaks = colorbreaksexp
colorbreaklabelsat = c(-2,-1,0,1,2,3)
colorbreaklabels = 10^colorbreaklabelsat
#cols = heat.colors(length(colorbreaks)) # for ppt
cols = gray(seq(from=0.1, to=1, length.out = length(colorbreaks)))
lwdlines = 3 # line weight on the graph: ppt
#axislwd=1 # line weight for axis and tickmarks: paper
axislwd =2# ppt
#fontsize = 12 # base font size: paper
fontsize = 22 # ppt
axiscex = 1 # expansion size for axis font: paper
axiscex = 0.75 # expansion size for axis font: ppt
labelcex = 1 # expansion size for label font


#title = paste("Genetic diversity (theta) \n(Philippines ",  km, "k, mu=", mu, ", stepwise model)",sep="") # main title for notebook
title = "Genetic diversity (theta)"

trellis.par.set(par.xlab.text=list(cex=labelcex),par.ylab.text=list(cex=labelcex),par.zlab.text=list(cex=labelcex), fontsize=list(text=fontsize))

contourplot(logthetaS ~ sigma+dens, data=he, region=TRUE, contour=FALSE, at=colorbreaks, col.regions = rev(cols), colorkey = list(labels=list(labels=colorbreaklabels, at = colorbreaklabelsat)), scales=list(log = 10, cex=axiscex), ylab="Effective density (De)", xlab = "Dispersal spread (sigma)", main = title, par.settings = list(axis.line = list(lwd = axislwd)), aspect="iso", 
	panel = function(at,contour,region,lwd,labels,...){ 
	panel.contourplot(contour=FALSE,at=colorbreaks,...) # plot the underlying colors
	panel.contourplot(at=c(-5,thetalog), labels = list(labels=c("",round(theta,2))), region=FALSE, contour=TRUE,lwd=lwdlines,label.style="align",...) # add a contour`
	panel.segments(log10(sigmamax),min(densexp),log10(sigmamax),Decritexp, lty=2, lwd=lwdlines)
	panel.segments(min(sigmaexp), Decritexp,log10(sigmamax), Decritexp, lty=2, lwd=lwdlines)
})
print(paste("Decrit =", Decrit))

# for ppt, without lines
contourplot(logthetaS ~ sigma+dens, data=he, region=TRUE, contour=FALSE, at=colorbreaks, col.regions = rev(heat.colors(length(colorbreaks))), colorkey = list(labels=list(labels=colorbreaklabels, at = colorbreaklabelsat)), scales=list(log = 10, cex=axiscex), ylab="Effective density (De)", xlab = "Dispersal spread (sigma)", main = title, par.settings = list(axis.line = list(lwd = axislwd)), aspect="iso", 
	panel = function(at,contour,region,lwd,labels,...){ 
	panel.contourplot(contour=FALSE,at=colorbreaks,...) # plot the underlying colors
})


#######################################################
## Plot De as function of sigma for a given theta and various mus

sigma = c(.1, 1, 5, 10, 50, 100, 500, 750, 1000, 2000, 3000, 5000, 10000)

k = 24000 # kms of reef in Phils, from Reefbase GIS
#mu = c(1e-4, 5e-4, 1e-3) # possible mutation rates
#mu = c(1e-6, 1e-4, 5e-4, 1e-3, 1e-2) # possible mutation rates
mu = c(5e-4) # possible mutation rates
#theta = 3 # observed theta
theta = 3.55 # observed theta 5/31/09 (missing 24,25,27)
a = 3.545 # Gaussian dispersal

sigmamax = 1000 # max sigma to consider
ymin = 0.001

if(!(sigmamax %in% sigma)){
	print("You didn't include sigmamax in sigma!!")
}

quartz(width=6, height=6)
#col = heat.colors(length(sigma))
col = "dark blue"
par(cex=1.3)
plot(.00001,.00001, ylim=c(ymin, 300), xlim=c(0.9*min(sigma), 1.1*max(sigma)), log="xy", ylab = "Effective density, De (individuals/km)", xlab = "Dispersal spread, sigma (km)", main=paste("Effective density in continous populations\n(theta = ", theta, ", k = ", k, ", a = ", a, ")", sep=""))
for(i in 1:length(mu)){
	De = theta/(4*mu[i]*sqrt(a*sigma*k))
	lines(sigma, De, col=col[i], lwd=2)
	Decrit = De[sigma==sigmamax]
	lines(x = c(sigmamax, sigmamax), y = c(ymin, Decrit), lty = 2, col="black")
	lines(x=c(0.9*min(sigma), 1000), y = c(Decrit, Decrit), lty=2)
	print(paste("mu = ", mu[i], ", Decrit = ", Decrit, sep=""))
}

legend("topright", legend = mu, fill=col, title="Mutation rate")




#######################################################
## Plot De as function of sigma for a given He and various mus

sigma = c(.1, 1, 5, 10, 50, 100, 500, 750, 1000, 2000, 3000, 5000, 10000)

k = 24000 # kms of reef in Phils, from Reefbase GIS
#mu = c(1e-4, 5e-4, 1e-3) # possible mutation rates
#mu = c(1e-6, 1e-4, 5e-4, 1e-3, 1e-2) # possible mutation rates
mu = c(5e-4) # possible mutation rates
He = 0.615 # observed heterozygosity
a = 3.545 # Gaussian dispersal

sigmamax = 1000 # max sigma to consider
ymin = 0.001

if(!(sigmamax %in% sigma)){
	print("You didn't include sigmamax in sigma!!")
}

quartz(width=6, height=6)
#col = heat.colors(length(sigma))
col = "dark blue"
par(cex=1.3)
plot(.00001,.00001, ylim=c(ymin, 300), xlim=c(0.9*min(sigma), 1.1*max(sigma)), log="xy", ylab = "Effective density, De (individuals/km)", xlab = "Dispersal spread, sigma (km)", main=paste("Effective density in continous populations\n(He = ", He, ", k = ", k, ", a = ", a, ")", sep=""))
for(i in 1:length(mu)){ # for each mutation rate
	De = 1/((1-He)^2*8*mu[i]*sqrt(a*sigma*k))-1/(8*mu[i]*sqrt(a*sigma*k))
	lines(sigma, De, col=col[i], lwd=2)
	Decrit = De[sigma==sigmamax]
	lines(x = c(sigmamax, sigmamax), y = c(ymin, Decrit), lty = 2, col="black")
	lines(x=c(0.9*min(sigma), 1000), y = c(Decrit, Decrit), lty=2)
	print(paste("mu = ", mu[i], ", Decrit = ", Decrit, sep=""))
}

legend("topright", legend = mu, fill=col, title="Mutation rate")



##################################################
### Plot LAMARC results
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090507")
lam = read.csv("LAMARC/Lamarc Summary 090531.csv", skip=2, header=T)

boxplot(lam$MPE)


dens = read.csv("../../../Surveys/surveys2009-03-31.density.csv")
rand = dens[dens$RandomSite==1,c("SiteNum", "densAPCL")]

lamdens = merge(lam, rand, by.x="Focal.pop", by.y = "SiteNum")
lamdens$region = "Cebu"
i = lamdens$Focal.pop == 19 | lamdens$Focal.pop == 20 | lamdens$Focal.pop == 22 | lamdens$Focal.pop == 23 | lamdens$Focal.pop == 24 | lamdens$Focal.pop == 25 | lamdens$Focal.pop == 27 
lamdens$region[i] = "Leyte"
lamdens$region = as.factor(lamdens$region)
plot(lamdens$densAPCL, lamdens$MPE, col=as.numeric(lamdens$region), pch=16)
plot(lamdens$densAPCL[lamdens$region=="Cebu"], lamdens$MPE[lamdens$region=="Cebu"], col=as.numeric(lamdens$region), pch=16)


plot(lamdens$Focal.pop, lamdens$MPE/max(lamdens$MPE), pch=16, type="b", ylim=c(0,1))
lines(lamdens$Focal.pop, lamdens$densAPCL/max(lamdens$densAPCL), pch=16, type="b", col="red")

summary(lm(lamdens$MPE ~ lamdens$densAPCL, subset=lamdens$region == "Cebu")) # Theta not related to densAPCL in Cebu (p = 0.11)
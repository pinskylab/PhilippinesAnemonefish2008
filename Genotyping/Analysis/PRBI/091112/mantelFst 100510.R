library(vegan)
library(smatr)

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/PRBI/091112/")

## Read in data

fst = read.csv("Genepop/PRBI_2009_11_12_fsts_E-W xsect.csv", row.names=1) # Genepop with all pops
geo = read.csv("PRBI_2010-05-10_geo.csv", row.names=1) 
geoalt2 = read.csv("PRBI_2010-05-10_counter_geo.csv", row.names=1) # measuring from 19 south of Bohol

## The APCL data for comparison
apclfst = read.csv("../../APCL/090507/genepop/Aclarkii_091028_all_1pop_fst.csv", header=TRUE, row.names=1)
apclgeo = read.csv("../../APCL/090507/Aclarkii_2009-05-14 geo.csv", row.names=1)


## Set up PRBI data
# linearize fst, remove diagonals, trim geo to those in fst
geo[upper.tri(geo, diag=T)] = NA
geoalt2[upper.tri(geoalt2, diag=T)] = NA
fst[upper.tri(fst, diag=T)] = NA
fstlin = as.matrix(fst/(1-fst))
	i = names(geo) %in% names(fst)
	geo = as.matrix(geo[i, i]) # trim out populations that we dont' use
	i = names(geoalt2) %in% names(fst)
	geoalt2 = as.matrix(geoalt2[i, i])

## Set up APCL data: linearize FSTs and split into Cebu and Leyte
apclgeo[upper.tri(apclgeo, diag=T)] = NA
apclfst[upper.tri(apclfst, diag=T)] = NA
apclfstlin = apclfst/(1-apclfst)

	cebupops = c('X7', 'X8', 'X9', 'X10', 'X11', 'X13', 'X14', 'X15', 'X16', 'X17')
	leytepops = c('X18', 'X19', 'X20', 'X22', 'X23', 'X24', 'X25', 'X26', 'X27')
	i = names(apclfstlin) %in% cebupops
	apclcebu_fst = as.matrix(apclfstlin[i,i])
	i = names(apclgeo) %in% cebupops
	apclcebu_geo = as.matrix(apclgeo[i,i])
	i = names(apclfstlin) %in% leytepops
	apclleyte_fst = as.matrix(apclfstlin[i,i])
	i = names(apclgeo) %in% leytepops
	apclleyte_geo = as.matrix(apclgeo[i, i])

##############
## Analysis ##
##############

# PRBI (all)
	plot(c(geo), c(fstlin), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West")
	l = line.cis(y=c(fstlin), x= c(geo))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1) # RMA linear model estimate

	mantel(geo, fstlin) # r = 0.117, p = 0.306


# PRBI w/out 19
	i = !(colnames(fstlin) %in% 'X19')
	j = !(colnames(geo) %in% 'X19')
	plot(c(geo[j,j]), c(fstlin[i,i]), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West w/out 19")
	l = line.cis(y=c(fstlin[i,i]), x= c(geo[j,j]))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1) # RMA linear model estimate

	# abline(a=l$coef[1], b=1/(4*52.7*4^2), col="grey", lwd=2) # slope expected with sigma=4km

	mantel(geo[j,j], fstlin[i,i]) # r = 0.533, p = 0.012

# PRBI Cebu only
	i = !(colnames(fstlin) %in% c('X1', 'X2', 'X19'))
	j = !(colnames(geo) %in% c('X1', 'X2', 'X19'))
	plot(c(geo[j,j]), c(fstlin[i,i]), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West w/out 19")
	l = line.cis(y=c(fstlin[i,i]), x= c(geo[j,j]))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1) # RMA linear model estimate

	mantel(geo[j,j], fstlin[i,i]) # r = 0.214, p = 0.285


# PRBI only 19, 7, 8, 9, 10, or 11 (genetic distance from 19 to 7 and 8 is bizarrely small)
	#pdf(paste('Figures/PRBI IBD by pop ', Sys.Date(), '.pdf', sep=""), width=8, height=5)
	xlims = c(0, 255)
	ylims = c(min(c(fstlin), na.rm=T), max(c(fstlin), na.rm=T))
	pchs=c(16, 13, 5, 2, 8, 1)
	cols=c('black', 'red', 'green', 'blue', 'orange', 'purple')
	ltys=2
	i = colnames(fstlin) %in% 'X19'
	j = colnames(geo) %in% 'X19'
	plot(c(geo[j,]), c(fstlin[i,]), pch=pchs[1], xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West by population", ylim=ylims, 	xlim=xlims, col=cols[1])
	l = line.cis(y=c(fstlin[i,]), x= c(geo[j,]))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1, lty=ltys) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X7'
	j = colnames(geo) %in% 'X7'
	points(c(geo[j,], geo[,j]), c(fstlin[i,], fstlin[,i]), pch=pchs[2], col=cols[2])
	l = line.cis(y=c(fstlin[i,], fstlin[,i]), x= c(geo[j,], geo[,j]))
	abline(a=l$coef[1], b=l$coef[2], col="red", lwd=1, lty=ltys) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X8'
	j = colnames(geo) %in% 'X8'
	points(c(geo[j,], geo[,j]), c(fstlin[i,], fstlin[,i]), pch=pchs[3], col=cols[3])
	l = line.cis(y=c(fstlin[i,], fstlin[,i]), x= c(geo[j,], geo[,j]))
	abline(a=l$coef[1], b=l$coef[2], col='green', lwd=1, lty=ltys) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X9'
	j = colnames(geo) %in% 'X9'
	points(c(geo[j,], geo[,j]), c(fstlin[i,], fstlin[,i]), pch=pchs[4], col=cols[4])
	l = line.cis(y=c(fstlin[i,], fstlin[,i]), x= c(geo[j,], geo[,j]))
	abline(a=l$coef[1], b=l$coef[2], col='blue', lwd=1, lty=ltys) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X10'
	j = colnames(geo) %in% 'X10'
	points(c(geo[j,], geo[,j]), c(fstlin[i,], fstlin[,i]), pch=pchs[5], col=cols[5])
	l = line.cis(y=c(fstlin[i,], fstlin[,i]), x= c(geo[j,], geo[,j]))
	abline(a=l$coef[1], b=l$coef[2], col='orange', lwd=1, lty=ltys) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X11'
	j = colnames(geo) %in% 'X11'
	points(c(geo[j,], geo[,j]), c(fstlin[i,], fstlin[,i]), pch=pchs[6], col=cols[6])
	l = line.cis(y=c(fstlin[i,], fstlin[,i]), x= c(geo[j,], geo[,j]))
	abline(a=l$coef[1], b=l$coef[2], col='purple', lwd=1, lty=ltys) # RMA linear model estimate

	legend('topright', bty='n', legend=c('19', '7', '8', '9', '10', '11'), pch=pchs, col=cols, lwd=1, lty=ltys)


## to output next 2 graphs to pdf:
	pdf(paste("Figures/PRBI IBD graphs ", Sys.Date(), ".pdf", sep=""), width=8, height=11)
	par(mfrow=c(2,1))
		# then run graph code, then run: dev.off()


# PRBI w/ and w/out 19
	#pdf(paste("Figures/PRBI IBD w and wout 19 ", Sys.Date(), ".pdf", sep=""), width=8, height=6)
	xlims = c(0, 255)
	ylims = c(min(c(fstlin), na.rm=T), max(c(fstlin), na.rm=T))
	i = !(colnames(fstlin) %in% 'X19')
	j = !(colnames(geo) %in% 'X19')
	plot(c(geo[j,j]), c(fstlin[i,i]), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West", ylim=ylims, xlim=xlims)
	l = line.cis(y=c(fstlin[i,i]), x= c(geo[j,j]))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1, lty=2) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X19'
	j = colnames(geo) %in% 'X19'
	points(c(geo[j,]), c(fstlin[i,]), pch=16, col='red', cex=1.2)
	l = line.cis(y=c(fstlin), x= c(geo))
	#abline(a=l$coef[1], b=l$coef[2], col='red', lwd=1, lty=1) # RMA linear model estimate

	legend('topleft', bty='n', legend=c('PRBI w/out #19', 'PRBI #19'), col=c('black', 'red'), pch=16)
	#dev.off()

# PRBI and APCL together (no Pop 19 PRBI)
	quartz(width=8, height=6)
	#pdf(paste("Figures/APCL IBD graphs ", Sys.Date(), ".pdf", sep=""), width=8, height=6)
	#pdf(paste("Figures/PRBI&APCL IBD graphs ", Sys.Date(), ".pdf", sep=""), width=8, height=6)
	xlims = c(0, 255)
	ylims = c(min(c(c(fstlin), c(apclcebu_fst), c(apclleyte_fst)), na.rm=T), max(c(c(fstlin), c(apclcebu_fst), c(apclleyte_fst)), na.rm=T))

	plot(c(apclcebu_geo), c(apclcebu_fst), pch=16, col='dark grey', cex=0.7, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="APCL vs. PRBI East-West w/out PRBI #19", ylim=ylims, xlim=xlims)
	l = line.cis(y=c(apclcebu_fst), x=c(apclcebu_geo))
	abline(a=l$coef[1], b=l$coef[2], col="dark grey", lwd=1, lty=1)

	points(c(apclleyte_geo), c(apclleyte_fst), pch=1, col='dark grey', cex=0.7)
	l = line.cis(y=c(apclleyte_fst), x=c(apclleyte_geo))
	abline(a=l$coef[1], b=l$coef[2], col="dark grey", lwd=1, lty=2)

	i = !(colnames(fstlin) %in% 'X19')
	j = !(colnames(geo) %in% 'X19')
	points(c(geo[j,j]), c(fstlin[i,i]), pch=16, col='red')
	l = line.cis(y=c(fstlin[i,i]), x= c(geo[j,j]))
	abline(a=l$coef[1], b=l$coef[2], col="red", lwd=1, lty=1) # RMA linear model estimate

	legend('topleft', legend=c('APCL Cebu', 'APCL Leyte', 'PRBI w/out #19'), pch = c(16,1,16), col=c('dark grey', 'dark grey', 'red'), pt.cex=c(0.7, 0.7,1), lty=c(1,2,1), bty='n')
	
	#legend('topleft', legend=c('APCL Cebu', 'APCL Leyte'), pch = c(16,1), col=c('dark grey', 'dark grey'), pt.cex=c(0.7, 0.7), lty=c(1,2), bty='n')

	# abline(a=l$coef[1], b=1/(4*52.7*4^2), col="grey", lwd=2) # PRBI slope expected with sigma=4km


	dev.off()


# PRBI and APCL together (Pop 19 PRBI colored separately)
	xlims = c(0, 255)
	ylims = c(min(c(c(fstlin), c(apclcebu_fst), c(apclleyte_fst)), na.rm=T), max(c(c(fstlin), c(apclcebu_fst), c(apclleyte_fst)), na.rm=T))
	i = !(colnames(fstlin) %in% 'X19')
	j = !(colnames(geo) %in% 'X19')
	plot(c(geo[j,j]), c(fstlin[i,i]), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="APCL vs. PRBI East-West", ylim=ylims, xlim=xlims)
	l = line.cis(y=c(fstlin[i,i]), x= c(geo[j,j]))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1, lty=1) # RMA linear model estimate

	i = colnames(fstlin) %in% 'X19'
	j = colnames(geo) %in% 'X19'
	points(c(geo[j,]), c(fstlin[i,]), pch=16, col='red')
	l = line.cis(y=c(fstlin), x= c(geo))
	abline(a=l$coef[1], b=l$coef[2], col='red', lwd=1, lty=1) # RMA linear model estimate

	points(c(apclcebu_geo), c(apclcebu_fst), pch=16, col='green', cex=0.5)
	l = line.cis(y=c(apclcebu_fst), x=c(apclcebu_geo))
	abline(a=l$coef[1], b=l$coef[2], col="green", lwd=1, lty=2)

	points(c(apclleyte_geo), c(apclleyte_fst), pch=16, col='blue', cex=0.5)
	l = line.cis(y=c(apclleyte_fst), x=c(apclleyte_geo))
	abline(a=l$coef[1], b=l$coef[2], col="blue", lwd=1, lty=2)

	legend('topleft', legend=c('PRBI', 'PRBI #19', 'APCL Cebu', 'APCL Leyte'), pch = 16, col=c('black', 'red', 'green', 'blue'), pt.cex=c(1,1,0.5, 0.5), lty=c(1,1,2,2), bty='n')




# Compare PRBI w/ two alternative distance matrices

	xlims = c(0, max(c(geo), c(geoalt2), na.rm=T))
	plot(c(geo), c(fstlin), pch=16, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West", xlim=xlims)
	l = line.cis(y=c(fstlin), x= c(geo))
	abline(a=l$coef[1], b=l$coef[2], col="black", lwd=1) # RMA linear model estimate

	points(c(geoalt2), c(fstlin), pch=2, xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI East-West", col='red', cex=1)
	l = line.cis(y=c(fstlin), x= c(geoalt2))
	abline(a=l$coef[1], b=l$coef[2], col="red", lwd=1) # RMA linear model estimate


	mantel(geo, fstlin) # r = 0.117, p = 0.306
	mantel(geoalt2, fstlin) # r = 0.4925, p = 0.016

## PRBI for PPT: white on black, w/out 19
	quartz(width=7, height=7)
	# pdf(width=7, height=7, file=paste('Figures/PRIB IBD wout19 PPT ', Sys.Date(), '.pdf', sep=''))
	col = "white"
	col.dot = 'cornflowerblue'
	par(cex=2, mar = c(4,3,3,0.5), mgp=c(2,0.5,0), omi=c(0,0,0,0), cex.axis=0.8, cex.lab=0.8, bg='black', bty='n', col.axis=col, col.lab=col, col=col, col.sub=col, col.main=col)
	i = !(colnames(fstlin) %in% 'X19')
	j = !(colnames(geo) %in% 'X19')
	xlims = range(geo[j,j], na.rm=T)
	ylims = range(fstlin[i,i], na.rm=T)
	
	plot(as.dist(geo[j,j]), as.dist(fstlin[i,i]), xlab="Geographic Distance (km)", ylab="Fst/(1-Fst)", main="PRBI", ylim=ylims, xlim = xlims, pch = 20, lwd=2, col=col.dot, xaxt='n', yaxt='n')
	axis(1, col=col, lwd=3)
	axis(2, col=col, lwd=3)
	
		# Regression line
		x = as.dist(geo[j,j])[1:length(as.dist(geo[j,j]))]
		y = as.dist(fstlin[i,i])[1:length(as.dist(fstlin[i,i]))]
		mod=sma(y ~ x)
		k = order(x)
		lines(x[k], mod$coef[[1]][1,1] + x[k]*mod$coef[[1]][2,1], lty=1, lwd=3, col=col.dot)

	dev.off()

	summary(mod)
	mantel(as.dist(fstlin[i,i]), as.dist(geo[j,j]), permutations=10000)

#################################################################
## Resampling approach to estimate 95% CI on dispersal distance
#################################################################

# From Rousset method: slope = 1/4Dsigma2
len = 10000 # how many iterations?

# RMA slopes and 95% CIs for calculation of log-normal SE (with and without pop 19)
slopes = data.frame(region=c("Wout19", "W19"), m = c(1.04068e-4, 8.39104835e-5), l95 = c(7.004779e-5, 5.67915025278822e-05), u95 = c(1.546115e-4, 0.0001239792738144356), sdlog = c(NA, NA))
for(i in 1:2){
	j = log(slopes$m[i])-log(slopes$l95[i])
	k = log(slopes$u95[i])-log(slopes$m[i])
	print(paste(j,k))
	slopes$sdlog[i] = mean(c(j,k)/1.96)
}
slopes


# Dataframe of means and sd for m, D, Ne/N, and De
	# Using RMA linear model regressions from R with approximate SE from log-normal CI using line.cis in smatr
priors = list(names=c("Wout19", "W19"), 
m = slopes$m, # IBD slopes
msdlog = slopes$sdlog, # IBD slope sdlog (for use in rlnorm())
D = c(52.7, 88.4), # census density
Dse = c(14.4, 37.8), # SE of D
nen = c(0.57, 0.57), # Ne/N ratio
nense = c(0.065, 0.065), # SE of nen
sigma = rep(NA, 14), # the point estimate
sigmas = vector("list", 14), # the sampled estimates
Des = vector("list", 14)) 

for(i in 1:length(priors$names)){
	D = rnorm(len, mean = priors$D[i], sd = priors$Dse[i])
	nen = rnorm(len, mean=priors$nen[i], sd = priors$nense[i])
	De = D*nen
	priors$Des[i] = list(De)
	priors$sigma[i] = sqrt(1/(4*priors$D[i]*priors$nen[i]*priors$m[i])) # the point estimate
	m = rlnorm(len, mean = log(priors$m[i]), sd = priors$msdlog[i])
	sigma = sqrt(1/(4*De*m))
	priors$sigmas[i] = list(sigma) 
}

# Print results
for(i in 1:length(priors$names)){
	temp = priors$sigmas[[i]]
	print(paste(priors$names[i], "  Point:", round(priors$sigma[i], digits=1), "  Mean:", round(mean(temp, na.rm=T), digits=1), "  Median:", round(median(temp, na.rm=T), digits=1), " Sd:", round(sd(temp, na.rm=T), digits=1), "  L95:", round(quantile(temp, 0.025, na.rm=T), digits=1), "  U95:", round(quantile(temp, 0.975, na.rm=T), digits=1)))
}



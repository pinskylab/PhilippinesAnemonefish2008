## Application of the Kaplan et al. 2006 Ecol Appl Dispersal Per Recruit (DPR) method


fsat = 0.35 #(p. 2252)



## Figure 2: single isolated reserve
#fneprfish = 0.2 # fishing leaves some behind
fneprfish = 0 # scorched-earth
mpa = 10 # width of reserve
shoulders = 100 # width of fished area on either side of reserve
fnepr = c(rep(fneprfish,shoulders),rep(1,mpa), rep(fneprfish,shoulders)) # Fraction of Natural Egg Per Recruit (reserve is center)

sat07 = eq4(0.7*mpa, fnepr)
sat09 = eq4(0.9*mpa, fnepr)
sat11 = eq4(1.1*mpa, fnepr)
sat16 = eq4(1.6*mpa, fnepr)
sat22 = eq4(2.2*mpa, fnepr)

# testing with eq6 for collapse
sat09_2 = eq6(0.9*mpa, fnepr, fsat, 100)
sat11_2 = eq6(1.1*mpa, fnepr, fsat, 100)
sat16_2 = eq6(1.6*mpa, fnepr, fsat, 100)
sat22_2 = eq6(2.2*mpa, fnepr, fsat, 100)


xlims=c(90,120)
ylims=c(0,1)
plot(0,0, xlim=xlims, ylim=ylims, xlab="Coastline", ylab="Fraction of natural settlement")
polygon(x=c(shoulders+1,shoulders+1,shoulders+mpa,shoulders+mpa), y=c(-1,2,2,-1), col="grey", border=NA) # the marine reserve
abline(h=fsat, col="black") # the critical level
cols = c("black", "red")
lines(sat07, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=2)
lines(sat09, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=2, col="green")
lines(sat11, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col="red")
lines(sat16, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=1)
lines(sat22, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=4)

# plot equilibrium level
lines(sat09_2, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col="green")
lines(sat11_2, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col="red")


## Figure 3a, b, and c: a system of reserves where they occupy 30% or 15% of coastline
fneprfish = 0 # fig 3a
fneprfish = 0.2 # fig 3b
mpa = 10
shoulders = 20 # fig 3a or 3b
shoulders = 57 # fig 3c
fnepr = c(rep(c(rep(fneprfish, shoulders), rep(1, mpa)), 4), rep(fneprfish, shoulders))

a07 = kaplan(0.7*mpa, fnepr, fsat, 100)
a09 = kaplan(0.9*mpa, fnepr, fsat, 100)
a11 = kaplan(1.1*mpa, fnepr, fsat, 100)
a16 = kaplan(1.6*mpa, fnepr, fsat, 100)
a22  = kaplan(2.2*mpa, fnepr, fsat, 100)


sat07 = eq4(0.7*mpa, fnepr)
sat09 = eq4(0.9*mpa, fnepr)
sat11 = eq4(1.1*mpa, fnepr)
sat16 = eq4(1.6*mpa, fnepr)
sat22 = eq4(2.2*mpa, fnepr)


xlims=c(shoulders*1.5+mpa,shoulders*3.5+mpa*3)
ylims=c(0,1)
plot(0,0, xlim=xlims, ylim=ylims, xlab="Coastline", ylab="Fraction of natural settlement")
polygon(x=c(shoulders*2+mpa+1,shoulders*2+mpa+1,shoulders*2+2*mpa,2*shoulders+2*mpa), y=c(-1,2,2,-1), col="grey", border=NA) # one marine reserve
polygon(x=c(shoulders*3+mpa*2+1,3*shoulders+2*mpa+1,3*shoulders+3*mpa,3*shoulders+3*mpa), y=c(-1,2,2,-1), col="grey", border=NA) # second marine reserve
abline(h=fsat, col="black") # the critical level
cols = c("red", "black") # black is sustainable, red is collapse
lines(sat07, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=2, col=cols[(sum(a07>0)>0)+1])
lines(sat09, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=2, col=cols[(sum(a09>0)>0)+1])
lines(sat11, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col=cols[(sum(a11>0)>0)+1])
lines(sat16, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=1, col=cols[(sum(a16>0)>0)+1])
lines(sat22, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=4, col=cols[(sum(a22>0)>0)+1])



## Figures for Phils 2008 Season Report
## Spatial scale: 10 "units" equal 1 km
## Narrow MPAs with varying FNEPR
fsat = 0.5
fneprfish1 = 0 
fneprfish2 = 0.5
fneprfish3 = 0.8 
fneprfish4 = 0.2
mpa = 10 # 1km mpa
nummpas = 20
shoulders = 90 # mpa every 10 km
#shoulders = 20 # mpa every 3 km
fnepr1 = c(rep(c(rep(fneprfish1, shoulders), rep(1, mpa)), nummpas), rep(fneprfish1, shoulders))
fnepr2 = c(rep(c(rep(fneprfish2, shoulders), rep(1, mpa)), nummpas), rep(fneprfish2, shoulders))
fnepr3 = c(rep(c(rep(fneprfish3, shoulders), rep(1, mpa)), nummpas), rep(fneprfish3, shoulders))
fnepr4 = c(rep(c(rep(fneprfish4, shoulders), rep(1, mpa)), nummpas), rep(fneprfish4, shoulders))

a1 = kaplan(11*mpa, fnepr1, fsat, 100)
a2 = kaplan(11*mpa, fnepr2, fsat, 100)
a3 = kaplan(11*mpa, fnepr3, fsat, 100)
a4 = kaplan(11*mpa, fnepr4, fsat, 100)

par(mfrow=c(1,2))
xlims=c(shoulders*1.5+mpa,shoulders*3.5+mpa*3)
xlims=c(0, length(fnepr1))
xlims = c(1000, 1500)
ylims=c(0,1)
plot(-1,-1, xlim=xlims, ylim=ylims, xlab="Coastline (km)", ylab="Fraction of natural settlement", xaxt="n")
axis(1, at=seq(1000, 1500, by=100), labels=seq(0, 50, by=10))
for(i in 1:nummpas){
	polygon(x=c(i*shoulders+mpa*(i-1)+1,i*shoulders+mpa*(i-1)+1,i*shoulders+i*mpa,i*shoulders+i*mpa), y=c(-1,2,2,-1), col="grey", border=NA) # one marine reserve
}
#abline(h=fsat, col="black") # the critical level

cols = c("red", "black") # black is sustainable, red is collapse
lines(a1, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=1, col=cols[(sum(a1>0)>0)+1])
lines(a2, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=2, col=cols[(sum(a2>0)>0)+1])
lines(a3, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col=cols[(sum(a3>0)>0)+1])
lines(a4+0.01, xlim=xlims, ylim=ylims, type="l", lwd=1, lty=2, col=cols[(sum(a4>0)>0)+1])
plot(-1, -1, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
legend("left", legend=c(0, 0.2, 0.5, 0.8, "Sanctuary"), lty=c(1, 2, 2, 1, 1), lwd=c(1, 1, 2, 2, 5), bty="n", title="FNEPR in\nfished areas", col=c(rep("black", 4), "grey"))


## Spatial scale: 10 "units" equal 1 km
## MPAs of various widths with FNEPR=0.2
fsat = 0.5
fneprfish5 = 0.2
mpa1 = c(20, 90, 120, 150, 150) # 2km, 9, 12, 15 mpa
nummpas = c(10, 10, 10, 10, 5)
shoulders1 = c(200,200,200,200, 500) # mpa spacing 20 km
fnepr5 = c(rep(c(rep(fneprfish5, shoulders1[1]), rep(1, mpa1[1])), nummpas), rep(fneprfish5, shoulders1[1]))
fnepr6 = c(rep(c(rep(fneprfish5, shoulders1[2]), rep(1, mpa1[2])), nummpas), rep(fneprfish5, shoulders1[2]))
fnepr7 = c(rep(c(rep(fneprfish5, shoulders1[3]), rep(1, mpa1[3])), nummpas), rep(fneprfish5, shoulders1[3]))
fnepr8 = c(rep(c(rep(fneprfish5, shoulders1[4]), rep(1, mpa1[4])), nummpas), rep(fneprfish5, shoulders1[4]))
fnepr9 = c(rep(c(rep(fneprfish5, shoulders1[5]), rep(1, mpa1[5])), nummpas[5]), rep(fneprfish5, shoulders1[5]))
fnepr5to8 = list(fnepr5, fnepr6, fnepr7, fnepr8, fnepr9)

a5 = kaplan(110, fnepr5, fsat, 100)
a6 = kaplan(110, fnepr6, fsat, 100)
a7 = kaplan(110, fnepr7, fsat, 100)
a8 = kaplan(110, fnepr8, fsat, 100)
a9 = kaplan(110, fnepr9, fsat, 100)
a5to8 = list(a5, a6, a7, a8, a9)

par(mfrow=c(4,1), oma=c(0, 0, 0, 0), mar=c(3, 4, 1, 1))
ylims=c(-0.1,1)
cols = c("red", "black") # black is sustainable, red is collapse
for(j in c(1,2,3,5)){
	len = length(fnepr5to8[[j]])
	x1 = round(5*len/12)
	x2 = round(7*len/12)
	xlims=c(x1, x2)
	#xlims=c(0,length(fnepr8))
	plot(-1,-1, xlim=xlims, ylim=ylims, xlab="", ylab="FNS", xaxt="n")
	axis(1, at=seq(x1, x2, by=100), labels=round(seq(0, (x2-x1)/10, by=10)))
	for(i in 1:nummpas[j]){
		polygon(x=c(i*shoulders1[j]+mpa1[j]*(i-1)+1,i*shoulders1[j]+mpa1[j]*(i-1)+1,i*shoulders1[j]+i*mpa1[j],i*shoulders1[j]+i*mpa1[j]), y=c(-1,2,2,-1), col="grey", border=NA) # one marine reserve
	}
	#abline(h=fsat, col="black") # the critical level

	lines(a5to8[[j]], xlim=xlims, ylim=ylims, type="l", lwd=1, lty=1, col=cols[(sum(a5to8[[j]]>0)>0)+1])
}

plot(-1, -1, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
legend("left", legend=c(0, 0.2, 0.5, 0.8, "MPA"), lty=c(1, 2, 2, 1, 1), lwd=c(1, 1, 2, 2, 5), bty="n", title="FNEPR in\nfished areas", col=c(rep("black", 4), "grey"))





## Spatial scale: 2 "units" equal 1 km
fsat = 0.5
fneprfish2 = 0.5
mpa = 2 # 1km mpa
nummpas = 10
shoulders2 = 98 # mpa every 50 km
fnepr2.2 = c(rep(c(rep(fneprfish2, shoulders2), rep(1, mpa)), nummpas), rep(fneprfish2, shoulders2))

a2.2 = kaplan(11*mpa, fnepr2.2, fsat, 100)

xlims=c(shoulders*1.5+mpa,shoulders*3.5+mpa*3)
xlims=c(0, length(fnepr2.2))
ylims=c(0,1)
plot(-1,-1, xlim=xlims, ylim=ylims, xlab="Coastline", ylab="Fraction of natural settlement")
for(i in 1:nummpas){
	polygon(x=c(i*shoulders+mpa*(i-1)+1,i*shoulders+mpa*(i-1)+1,i*shoulders+i*mpa,i*shoulders+i*mpa), y=c(-1,2,2,-1), col="grey", border=NA) # one marine reserve
}
abline(h=fsat, col="black") # the critical level

cols = c("red", "black") # black is sustainable, red is collapse
lines(a2.2, xlim=xlims, ylim=ylims, type="l", lwd=2, lty=1, col=cols[(sum(a2.2>0)>0)+1])





#### Functions

laplace <- function(x,y,a){
	p = exp(-abs(x-y)/a)/(2*a)
	return(p)	
}

# test of habitat saturation condition (Eq. 4)
# a: mean dispersal distance in a laplacian kernel
# fnepr: fraction of natural eggs-per-recruit at each point along the habitat
# fsat: fraction of natural settlement needed to saturate post-recruitment habitat (0.35?)
eq4 = function(a, fnepr){
	len = length(fnepr)
	sat = numeric(len) 
	for(x in 1:len){ # calc recruitment at each point
		sum = 0
		for(y in (1:len)){ # calc recruitment from each point to here
			sum = sum + laplace(x,y,a)*fnepr[y]
		}
		sat[x] = sum
	}
	return(sat)
}

# Recruits relationship, given num settlers (Fig. 1)
# A hockey-stick relationship
# fns: fraction of natural settlement (can be a vector)
# fsat: fraction of natural settlement needed to saturate the habitat
#
# Fraction of natural recruitment = 1 for fns > fsat
# FNR declines linearly to 0 for fns = 0   
recruits = function(fsat, fns){
	i = which(fns >= fsat)
	j = which(fns < fsat)
	fnr = numeric(length(fns))
	if(length(i)>0) fnr[i] = 1
	if(length(j)>0) fnr[j] = 1/fsat * fns[j]
	return(fnr)
}

# DPRn: an iterative approach to determining the equilibrium population
# Eq. 6
# maxiter: maximum number of iterations to do (default 100)
# stops if fnr < fsat everywhere and sets fns to 0
eq6 = function(a, fnepr, fsat, maxiter=100){
	fns_old = eq4(a, fnepr) # fraction of natural settlement
	fns_new = eq4(a, recruits(fsat, fns_old)*fnepr) # fraction of natural settlement
	i = 0
	while((sum(fns_new != fns_old)>0) & (sum(fns_new > fsat)>0) & (i < maxiter)){
		fns_old = fns_new
		fns_new = eq4(a, recruits(fsat, fns_old)*fnepr) # fraction of natural settlement
		i = i+1
		if(i %% 10 == 0) print(i)
	}
	if(sum(fns_new > fsat)==0) fns_new = rep(0, length(fnepr))
	if(i == maxiter) warning("eq6: reached maxiter, may need to increase it")
	
	return(fns_new)
}

# Combination of Eq. 4 and Eq. 6: moves to 6 if fns < fsat at some locations where fnepr > 0
kaplan = function(a, fnepr, fsat, maxiter=100){
	fns = eq4(a, fnepr)
	ind = fnepr > 0 # places where fnepr is non-zero (places we care about)
	if(all(fns[ind] > fsat)){
		return(fns)	
	} else {
		fns = eq6(a, fnepr, fsat, maxiter)	
		return(fns)
	}
}
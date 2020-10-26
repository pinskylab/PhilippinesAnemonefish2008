# Jackknife the Arlequin heterozygosity estimates
# from http://www.math.wustl.edu/~sawyer/handouts/jackknife.pdf


setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/090507/Aclarkii_2009-05-08 CebuLeyte onepop.res/")

he = read.csv("Heterozygosity.csv")
n = length(he$Locus)

ps = numeric(n) # pseudovalues
Xn = n*mean(he$He)
for(i in 1:length(jack)){
	ps[i] = Xn - (n-1)*mean(he$He[-i])	
}

sum(ps)/n
sd(ps)/sqrt(n)

# or the basic way (same answer)
mean(he$He)
sd(he$He)/sqrt(n)

# produces strange results... jackknife the theta estimate?
he$theta = 1/(2*(1-he$He)^2)-1/2
mean(he$theta)
sd(he$theta)/sqrt(n)

# give this function a list of He, will return theta
returntheta = function(h){
	hbar = mean(h)
	return(1/(2*(1-hbar)^2)-1/2)
}

ps = numeric(n) # pseudovalues
Xn = n*returntheta(he$He)
for(i in 1:length(jack)){
	ps[i] = Xn - (n-1)*returntheta(he$He[-i])	
}

psbar = sum(ps)/n
psvar = 1/(n-1)*sum((ps - psbar)^2)
sqrt(psvar)/sqrt(n)

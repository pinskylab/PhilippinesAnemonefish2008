# plot data from MLNE Summary 090513.xls
require(car)

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/APCL/090507/")

mne = read.csv("MNE/Aclarkii_2009-05-13_MLNEfocalTopTwo/MLNE Summary 090513.csv") # MNe-All method
#mne = read.csv("MNE/Aclarkii_2009-06-29_MLNEfocalTopTwoSurrTwo/MLNE-Flanking 090629.csv") # MNe-Flanking method

mne$nlog10 = log10(mne$ML_Ne)
mne$region = NA
mne$region[mne$Pop<=17] = "Cebu"
mne$region[mne$Pop>=18] = "Leyte"

# summary stats
aggregate(mne$ML_Ne, by=list(mne$region), mean)$x # mean adults
as.numeric(aggregate(mne$ML_Ne, by=list(mne$region), sd)$x/sqrt(table(mne$region))) # SE adults

aggregate(mne$ML_Ne, by=list(mne$region), median)$x # median adults
as.numeric(aggregate(mne$ML_Ne, by=list(mne$region), sd)$x/sqrt(table(mne$region))*1.253) # SE adults for the median (analytical)

aggregate(mne$ML_Ne, by=list(mne$region), median)$x/25 # median adults/km
as.numeric(aggregate(mne$ML_Ne, by=list(mne$region), sd)$x/sqrt(table(mne$region))*1.253/25) # SE adults for the median (analytical)

# bootstrap CI for the median
len = 10000
medscebu = numeric(len)
medsleyte = numeric(len)
for(i in 1:len){
	medscebu[i] = median(sample(mne$ML_Ne[mne$region=="Cebu"], sum(mne$region=="Cebu"), replace=TRUE))
	medsleyte[i] = median(sample(mne$ML_Ne[mne$region=="Leyte"], sum(mne$region=="Leyte"), replace=TRUE))
}
#hist(medscebu)
#hist(medsleyte)
quantile(medscebu, probs=c(0.025, 0.975)) # 95% CI for median adults (percentile method)
quantile(medsleyte, probs=c(0.025, 0.975))
quantile(medscebu, probs=c(0.025, 0.975))/25 # 95% CI for median adults/km (percentile method)
quantile(medsleyte, probs=c(0.025, 0.975))/25

	abs(median(mne$ML_Ne[mne$region=="Cebu"])-quantile(medscebu, probs=c(0.025, 0.975))) # distance to upper and lower CI: assymetrical!
	abs(median(mne$ML_Ne[mne$region=="Leyte"])-quantile(medsleyte, probs=c(0.025, 0.975)))

	log(abs(median(mne$ML_Ne[mne$region=="Cebu"])-quantile(medscebu, probs=c(0.025, 0.975)))) # distance to upper and lower CI: less assymetrical on the log scale (but still not great)
	log(abs(median(mne$ML_Ne[mne$region=="Leyte"])-quantile(medsleyte, probs=c(0.025, 0.975))))

	mean(log(abs(median(mne$ML_Ne[mne$region=="Cebu"])-quantile(medscebu, probs=c(0.025, 0.975))))) # mean log(distance) to upper and lower CI: for use in log-normal approximation to error
	mean(log(abs(median(mne$ML_Ne[mne$region=="Leyte"])-quantile(medsleyte, probs=c(0.025, 0.975)))))


	# 95% CI Lunneborg's method (see http://www.uvm.edu/~dhowell/StatPages/Resampling/BootstMeans/bootstrapping_means.html)
	median(mne$ML_Ne[mne$region=="Cebu"])+(median(mne$ML_Ne[mne$region=="Cebu"])-quantile(medscebu, probs=c(0.025, 0.975))) 
	median(mne$ML_Ne[mne$region=="Leyte"])+(median(mne$ML_Ne[mne$region=="Leyte"])-quantile(medsleyte, probs=c(0.025, 0.975)))



# bootstrap CI for the median using package boot
library(boot)

# simple approach
samplemedian <- function(x, d) { # have to set up a function that boot can use (uses indices d)
  return(median(x[d]))
}
cboot = boot(as.numeric(mne$ML_Ne[mne$region=="Cebu"]), statistic=samplemedian, R=10000, sim="ordinary")
lboot = boot(as.numeric(mne$ML_Ne[mne$region=="Leyte"]), statistic=samplemedian, R=10000, sim="ordinary")

plot(cboot)
plot(lboot)

boot.ci(cboot, conf=0.95)
boot.ci(lboot, conf=0.95)

boot.ci(cboot, conf=0.95, type="perc") # just the percentile CI
boot.ci(lboot, conf=0.95, type="perc")

	# or use studentized t-intervals, but needs bootstrap estimate of variance
	# see http://www.psychology.mcmaster.ca/bennett/boot09/percentileT.pdf
	boot.median <- function(theData,ids){
		boot.sample <- theData[ids] # bootstrapped sample
		m = median(boot.sample)
		n = length(boot.sample)
		R = 200 # number of bootstrap-within-a-bootstrap replicates to do
		# create R x n matrix of resampled values (a for loop would be too slow)
		reboot.samples = matrix(data=sample(boot.sample, size=R*n, replace=T), nrow=n)
		# calculate median of each n=20 column
		reboot.medians = apply(reboot.samples, MARGIN=2, median)
		# calculate variance of medians
		m.var = var(reboot.medians)
		results = c(m, m.var)
		return(results)
	}
	cboot2 = boot(as.numeric(mne$ML_Ne[mne$region=="Cebu"]), statistic=boot.median, R=1000, sim="ordinary")
	lboot2 = boot(as.numeric(mne$ML_Ne[mne$region=="Leyte"]), statistic=boot.median, R=1000, sim="ordinary")
	
	boot.ci(cboot2)
	boot.ci(lboot2)

#histograms
par(mfrow=c(1,2))
hist(mne$ML_Ne[mne$region=="Cebu"], breaks=20)
hist(mne$ML_Ne[mne$region=="Leyte"], breaks=20)

#histograms
par(mfrow=c(1,2))
hist(mne$nlog10[mne$region=="Cebu"], breaks=20)
hist(mne$nlog10[mne$region=="Leyte"], breaks=20)


# boxplots
quartz(width=8, height=8)
par(cex=1.7, cex.lab = 1.5, cex.axis=1)
boxplot(nlog10 ~ region, data=mne, range=0, yaxt="n", main = "Effective size", col="dark grey")
power.axis(power=0, base=10, at = c(1,10,50, 100,250, 500, 1000), axis.title="Ne", side=2, cex=1.5*1.3)

# for ppt (log10 axis)
quartz(width=8, height=8)
par(omi = c(0,0.1,0,0))
boxplot(nlog10 ~ region, data=mne, range=0, yaxt="n", main = "Effective size", col="white", boxwex=0.7, staplewex=0.2, outwex=0.5, lwd=3, cex.axis=2, cex.main=2)
power.axis(power=0, base=10, at = c(1,10,50, 100,250, 500, 1000), axis.title="Ne", side=2, cex=2)
abline(h = log10(203))
abline(h = log10(296))

par(omi = c(0,0,0,0), mgp=c(4,1,0), mar=c(4,7,4,2))
boxplot(n ~ region, data=x, range=0, main = "Effective size", ylab="Ne", col="dark grey", boxwex=0.7, staplewex=0.2, outwex=0.5, lwd=3, cex.axis=2.3, cex.main=3.2, cex.lab=3.2, log="y")

# normal axis
quartz(width=8, height=8)
par(omi = c(0,0.1,0,0))
boxplot(ML_Ne ~ region, data=mne, range=0, yaxt="a", main = "Effective size", col="white", boxwex=0.7, staplewex=0.2, outwex=0.5, lwd=3, cex.axis=2, cex.main=2)
power.axis(power=0, base=10, at = c(1,10,50, 100,250, 500, 1000), axis.title="Ne", side=2, cex=2)
abline(h = log10(203))
abline(h = log10(296))


t.test(n ~ region, data=x)

median(x$n)
mean(x$n)/25
sd(x$n)/length(x$n)/25


require(vioplot)
a = mne$ML_Ne[mne$region == "Cebu"]
b = mne$ML_Ne[mne$region == "Leyte"]

vioplot(a,b)
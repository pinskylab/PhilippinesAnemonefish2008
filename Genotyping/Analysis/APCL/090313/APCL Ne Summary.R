# plot data from APCL Ne Summary.xls
require(car)

# MNe2.0 for one-pop
a = c(87, 94, 267, 1008, 62, 63, 91, 90, 139, 125) # Cebu Ne from MNe2.0

b = c(1790, 286, 62, 95, 279, 52, 51, 106) # Leyte

x = data.frame(n = c(a,b), nlog10 = log10(c(a,b)),  region = c(rep("Cebu", length(a)), rep("Leyte", length(b))))

#quartz(width=3, height=3)
quartz(width=8, height=8)
#par(cex=.7, cex.lab = 1.4, cex.axis=0.95)
par(cex=1.7, cex.lab = 1.5, cex.axis=1)

boxplot(nlog10 ~ region, data=x, range=0, yaxt="n", main = "Effective size", col="dark grey")
power.axis(power=0, base=10, at = c(1,10,50, 100,250, 500, 1000), axis.title="Ne", side=2, cex=1.5*1.3)

# for ppt
par(omi = c(0,0.1,0,0))
boxplot(nlog10 ~ region, data=x, range=0, yaxt="n", main = "Effective size", col="dark grey", boxwex=0.7, staplewex=0.2, outwex=0.5, lwd=3, cex.axis=2, cex.main=2)
power.axis(power=0, base=10, at = c(1,10,50, 100,250, 500, 1000), axis.title="Ne", side=2, cex=2)

par(omi = c(0,0,0,0), mgp=c(4,1,0), mar=c(4,7,4,2))
boxplot(n ~ region, data=x, range=0, main = "Effective size", ylab="Ne", col="dark grey", boxwex=0.7, staplewex=0.2, outwex=0.5, lwd=3, cex.axis=2.3, cex.main=3.2, cex.lab=3.2, log="y")



t.test(n ~ region, data=x)

median(x$n)
mean(x$n)/25
sd(x$n)/length(x$n)/25
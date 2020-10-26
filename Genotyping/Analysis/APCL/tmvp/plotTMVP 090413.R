setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys/Genotype Analysis/tmvp/")
require(locfit)
source("../../../../ABC_distrib/loc2plot_d.r") # Beaumont's functions for 1D and 2D data

tmvp = read.table("25/Aclarkii_2009-04-09_TMVP25_100K,1M.out")
names(tmvp) = c("x", "logL", "N1", "N2")
dim(tmvp)

plot(tmvp$N1, tmvp$logL)
plot(tmvp$N1, tmvp$logL, xlim=c(0,10000))

loc1stats(tmvp$N1, 0.05)

linestats(tmvp$N1, tmvp$logL, 0.05)

fit = locfit(logL ~ lp(N1), data=tmvp)
summary(fit)
par(mfrow=c(1,2))
plot(fit, get.data=FALSE, band="global", ylim=c(-251.935, -251.915))
plot(fit, get.data=TRUE, band="global", ylim=c(-251.935, -251.915))

plot(density(tmvp$N1))

xx = seq(0,1000000, by=10)
yy = predict(fit, newdata=xx)
plot(xx,yy)
xx[which.min(yy)]
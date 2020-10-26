setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys/Genotype Analysis/tmvp/")
require(locfit)
source("../../../../ABC_distrib/loc2plot_d.r") # Beaumont's functions for 1D and 2D data

tmvp = read.table("25/Aclarkii_2009-04-09_TMVP25_20K,10k.out")
tmvp = read.table("25/Aclarkii_2009-04-09_TMVP25_100K,1M.out")
names(tmvp) = c("x", "logL", "N1", "N2")
dim(tmvp)


# Returns mode and 95% CI for N1
loc1stats(tmvp$N1, 0.05) # max of 10k or 1M still isn't high enough

par(mfrow=c(1,2))
plot(locfit(~tmvp$N1))
hist(tmvp$N1)


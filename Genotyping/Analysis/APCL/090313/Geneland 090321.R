setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping")

library(Geneland)

geno <- read.table("Genotypes/Aclarkii_2009-03-20_Gnlndgen.csv") 
coord <- read.table("Genotypes/Aclarkii_2009-03-20_Gnlndgeo.csv") 
coord[1:10,] 
#fix(coord) 
#fix(geno) 
dim(coord) 
dim(geno)

plot(coord,xlab="Eastings",ylab="Northings",asp=1)

# Multiple chains: Run MCMC and process
	nrun <- 10 
	burnin <- 200 
	for(irun in 1:nrun) { 
		path.mcmc <- paste("./Analysis/090313/Geneland/correlated/",irun,"/",sep="")
		dir.create(path.mcmc)
		MCMC(coordinates=coord, genotypes=geno, varnpop=TRUE, npopmax=10, spatial=TRUE, freq.model="Correlated", nit=100000, thinning=100, path.mcmc=path.mcmc) 

		## MCMC postprocessing 
		PostProcessChain(coordinates=coord, genotypes=geno, path.mcmc=path.mcmc, nxdom=100, nydom=100, burnin=burnin) 
	} 

	## Computing average posterior probability 
	## with a burnin of 200 (* 100) iterations 
	lpd <- rep(NA, nrun) 
	for(irun in 1:nrun) { 
		path.mcmc <- paste("./Analysis/090313/Geneland/correlated/",irun,"/",sep="") 
		path.lpd <- paste(path.mcmc,"log.posterior.density.txt",sep="") 
		lpd[irun] <- mean(scan(path.lpd)[-(1:burnin)]) 
	} 
	## Runs sorted by decreasing average posterior probability: 
	order(lpd,decreasing=TRUE) 


# Plots
Plotnpop(path.mcmc=mcmcpath, burnin=200) 
PosteriorMode(coordinates=coord, path.mcmc=mcmcpath)

for(irun in 1:nrun){
	quartz()
	path.mcmc <- paste("./Analysis/090313/Geneland/correlated/",irun,"/",sep="")
	Plotnpop(path.mcmc=path.mcmc, burnin=200) 
}

# Fstats
Fstat.output(genotypes=geno,path.mcmc=mcmcpath) 


# A long run
		path.mcmc <- paste("./Analysis/090313/Geneland/correlated/10mil/",sep="")
		#dir.create(path.mcmc)
		MCMC(coordinates=coord, genotypes=geno, varnpop=TRUE, npopmax=10, spatial=TRUE, freq.model="Correlated", nit=10000000, thinning=10000, path.mcmc=path.mcmc) 

		## MCMC postprocessing 
		PostProcessChain(coordinates=coord, genotypes=geno, path.mcmc=path.mcmc, nxdom=100, nydom=100, burnin=200) 
		
		Plotnpop(path.mcmc=path.mcmc, burnin=200) 

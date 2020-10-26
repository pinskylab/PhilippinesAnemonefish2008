###############################################
## Output to MLNE2: All samples (TopTwo on anemone >= 8cm vs. lower quartile, some pops as focal, all others as source)



NeMax = 10000

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	


genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"

	
	focal = 1

	# Write to file
	mlnefile = file(paste("MLNE/Aclarkii_", Sys.Date(), "_MLNEfocal", prefix, paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 1 (open pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(1, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(0, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# num alleles per locus (only for adults and juvs)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = genowide$AdJuv == "Ad" | genowide$AdJuv == "Juv"
	for(i in cols){ # for each locus
		j = union(genowide[a,i], genowide[a,i+1])
		j = j[!is.na(j)]
		alleles[match(i,cols)] = list(sort(j))
		outline = c(outline, length(j))
	}
	cat(paste(outline, collapse=","),file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
			
	# num samples from focal pop
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# generations when samples were taken
	cat("0,1", file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	
	sites = matchall(focal, genowide$SiteNum)

	# For adults in focal pops:
	# Numbers of copies of each allele, at each locus, from each sample of the focal population
	cols = grep("A1.", names(genowide))
	a = intersect(which(genowide$AdJuv == "Ad"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), sep=""))
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	# For juvs in focal pop
	cols = grep("A1.", names(genowide))
	a = intersect(which(genowide$AdJuv == "Juv"), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Juvs: ", length(a), sep=""))
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}

	# For source pop (use all samples not from focal)
	cols = grep("A1.", names(genowide))
	a = setdiff(1:length(genowide$SiteNum), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Source: ", length(a), sep=""))
	for(i in cols){ # for each locus
		thesealleles = sort(c(genowide[a,i], genowide[a,i+1]))
		allalleles = alleles[[match(i,cols)]]
		outline = numeric()
		for(j in 1:length(allalleles)){
			outline = c(outline, sum(thesealleles == allalleles[j]))
		}
		cat(paste(outline, collapse=","), file=mlnefile, sep=",", append=T)
		cat("\n", file= mlnefile, append=T)
	}
	
	# number of starting points
	cat("1", file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
		
	
	close(mlnefile)	

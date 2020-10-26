setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
genowide=read.csv("PRBI_genowide_2009-11-01.csv", row.names=1)

########################
# Format for Dropout

	a1cols = grep("A1.", names(genowide), fixed=TRUE)
	locusnames = gsub("A1.", "", names(genowide)[a1cols], fixed=T)
	numindivs = dim(genowide)[1]
	allelesep = "|"

	# Write to file
	file = file(paste("../Datafiles/Dropout/PRBI_dropout_", Sys.Date(), ".txt", sep=""))
	open(file, "w")

	# Header: num capture sessions, num loci, "!" for locus names, then locus names on each line
	outline = c(1, length(a1cols), "!")
	cat(outline, file=file, sep=" ", append=F)
	cat("\n", file=file, append=T)
	outline = paste(locusnames, collapse="\n")
	cat(outline, file=file, sep=" ", append=T)
	cat("\n", file=file, append=T)
	
	prefix = paste("/", genowide$Sample, " 1 ", sep="") # start each line with sample name and capture session (1)
	# collapse locus alleles
	out = cbind(prefix, paste(formatC(genowide[,a1cols[1]], width=3, flag="0"), formatC(genowide[,(a1cols[1]+1)], width=3, flag="0"), sep=allelesep))
	for(i in a1cols[2:length(a1cols)]){ # for remaining loci
		out = cbind(out, paste(formatC(genowide[,i], width=3, flag="0"), formatC(genowide[,(i+1)], width=3, flag="0"), sep=allelesep))
	}
	out[out==paste(" NA", allelesep, " NA", sep="")] = paste("000", allelesep, "000", sep="")

	for(i in 1:dim(out)[1]){
		cat(paste(out[i,], collapse=" "), file=file, append=T)
		cat("\n", file= file, append=T)
	}	

	close(file)	




###############################################
## Genepop: Output each pop (used by LDNE, Microchecker)

	a1cols = grep("A1.", names(genowide), fixed=TRUE)
	locusnames = gsub("A1.", "", names(genowide)[a1cols], fixed=T)
	numindivs = dim(genowide)[1]
	pops = sort(unique(genowide$SiteNum))

	# Write all pops to files
	#cols = grep("A1.", names(genowide))
	#loci = unlist(strsplit(names(genowide)[cols], "A1.", fixed=T))[seq(2,26, by=2)]

	# Write to file
	filename = paste("../Datafiles/Genepop/PRBI_genepop_", Sys.Date(), ".gen", sep="")
	file = file(filename)
	open(file, "w")

	# Header
	cat("Title line: PRBI all pops", file= file, sep=",", append=F)
	cat("\n", file= file, append=T)
	cat(paste(locusnames, collapse="\n"), file= file, append=T)
	cat("\n", file= file, append=T) 

	for(j in 1:length(pops)){
		k = which(genowide$SiteNum==pops[j])
		if(length(k)>0){
			cat("Pop", file= file, sep="", append=T)
			cat("\n", file= file, append=T)
		
			# collapse locus alleles with ""
			out = paste("     ", rep(pops[j], length(k)), ",", sep="") # name each indiv according to population name
			out = cbind(out, paste(formatC(genowide[k,a1cols[1]], width=3, flag="0"), formatC(genowide[k,(a1cols[1]+1)], width=3, flag="0"), sep=""))
			for(c in a1cols[2:length(a1cols)]){
				out = cbind(out, paste(formatC(genowide[k,c], width=3, flag="0"), formatC(genowide[k,(c+1)], width=3, flag="0"), sep=""))
			}
			out[out==" NA NA"] = "000000"
			for(c in 1:dim(out)[1]){
				cat(paste(out[c,], collapse=" "), file=file, append=T)
				cat("\n", file= file, append=T)
			}	
		}
	}
	close(file)	



###############################################
## Genepop: Output each individual w/ UTM X,Y (for IBD by individual analyses in genepop)

	a1cols = grep("A1.", names(genowide), fixed=TRUE)
	locusnames = gsub("A1.", "", names(genowide)[a1cols], fixed=T)
	numindivs = dim(genowide)[1]
	pops = sort(unique(genowide$SiteNum))

	#whichindivs = 1:numindivs; suff="" # to output all individuals
	whichindivs = which(genowide$Region=="Cebu"); suff="cebu_" # output only Cebu individuals

	# Write to file
	filename = paste("../Datafiles/Genepop/PRBI_genepopIBD_", suff, Sys.Date(), ".gen", sep="")
	file = file(filename)
	open(file, "w")

	# Header
	cat("Title line: PRBI all pops by indiv for IBD", file= file, sep=",", append=F)
	cat("\n", file= file, append=T)
	cat(paste(locusnames, collapse="\n"), file= file, append=T)
	cat("\n", file= file, append=T) 

	for(j in whichindivs){
		cat("Pop", file= file, sep="", append=T)
		cat("\n", file= file, append=T)
		
		# collapse locus alleles with ""
		out = paste(round(genowide$x_utm[j]), " ", round(genowide$y_utm[j]), ",", sep="") # name the indiv by its UTM coords
		out = cbind(out, paste(formatC(genowide[j,a1cols[1]], width=3, flag="0"), formatC(genowide[j,(a1cols[1]+1)], width=3, flag="0"), sep="")) # add genotypes
		for(c in a1cols[2:length(a1cols)]){
			out = cbind(out, paste(formatC(genowide[j,c], width=3, flag="0"), formatC(genowide[j,(c+1)], width=3, flag="0"), sep=""))
		}
		out[out==" NA NA"] = "000000"
		cat(paste(out, collapse=" "), file=file, append=T)
		cat("\n", file= file, append=T)
	
	}

	close(file)	


########################
# Output for create w/ UTM coords (conversion to spagedi)

	# Cebu only?
	k = which(genowide$Region=="Cebu")
	suff = "cebu_"
		
	# All
	k = 1:dim(genowide)[1]
	suff=""
	
	# Write to file
	file = file(paste("../Datafiles/CREATE/PRBI_create_", suff, Sys.Date(), ".csv", sep=""))
	open(file, "w")
	

	i = grep("A1", names(genowide))
	j = names(genowide)[i[1]:(i[length(i)]+1)]
	j = paste(j, collapse=",")
	j = gsub(".", "_", j, fixed=T)
	outline = paste("Indiv,X,Y,",j, sep="") 
	cat(outline, file=file, append=T)
	cat("\n", file=file, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[k,c(grep("Sample",widenames), grep("x_utm",widenames), grep("y_utm",widenames), grep("A[[:digit:]]",widenames))]
	write.table(outfile, file=file, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(file)	


########################
# Output for create w/ Adult-Juvenile (conversion to FaMoz and Cervus)

	# add Ad/juv info: Adults are top two on anemone if >=8cm, juvs are everything else
	cutoff = 8
	genowide$AdJuv = 0
	i = genowide$TopTwo & genowide$Size > cutoff
	genowide$AdJuv[i] = 0 # Adults
	genowide$AdJuv[!i] = 1 # Juvs

	# remove pops 18 and 7?
	#genowide = genowide[genowide$SiteNum != 7 & genowide$SiteNum != 18,]

	
	# Write to file
	file = file(paste("PRBI_create_cut", cutoff, "_", Sys.Date(), ".csv", sep=""))
	open(file, "w")
	

	i = grep("A1", names(genowide))
	j = names(genowide)[i[1]:(i[length(i)]+1)]
	j = paste(j, collapse=",")
	j = gsub(".", "_", j, fixed=T)
	outline = paste("Pop,Indiv,Cohort,",j,sep="") 
	cat(outline, file=file, append=T)
	cat("\n", file=file, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[,c(grep("SiteNum",widenames),grep("Sample",widenames),grep("AdJuv", widenames), grep("A[[:digit:]]",widenames))]
	write.table(outfile, file=file, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(file)	

	filename = paste("PRBI_genowide_", Sys.Date(), ".csv", sep="")
	# write.csv(genowide, filename)





###############################################
## Output to MLNE2: All samples (upper quartile vs. lower quartile, some pops as focal, all others as source)
source("../Analysis/matchall.R")	

# sets of 4
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27))
prefix="TopTwo"

# sets of 4, 3, 2, and 1 (for testing only)
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27), m = c(7,8,9), n = c(8,9,10), o=c(9,10,11), p=c(10,11,13), q=c(11,13,14), r=c(13,14,15), s=c(14,15,16), t=c(15,16,17), u=c(18,19,20), v=c(19,20,22), w=c(20,22,23), x=c(22,23,24), y=c(23,24,25), z=c(24,25,27), aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27), qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)

# sets of 1,2, and 3, all using pop 25 (determined from sets above)
focalpops = list(a = 25, b = c(24,25), c=c(23,24,25))
prefix="TopTwo"

# sets of 2
focalpops = list(aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27))
prefix="TopTwo"

# sets of 1
focalpops = list(oo = 1, pp = 2, qq=7, rr=8, ss=9, tt=10, uu=11, bbb=19)

NeMax = 10000

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75), na.rm=TRUE)

first = genowide$Size <= quants[1]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[fourth_inclusive] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

# How many breeding adults (upper quartile) and juvs (lowest quartile) in each population group?
for(rep in 1:length(focalpops)){	
	focal = focalpops[[rep]]	
	sites = matchall(focal, genowide$SiteNum)
	a = intersect(which(genowide$sizeclass == 4), sites)
	b = intersect(which(genowide$sizeclass == 1), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

# Write all groups to files
for(rep in 1:length(focalpops)){
	
	focal = focalpops[[rep]]

	# Write to file
	mlnefile = file(paste("../Datafiles/MLNE/PRBI_", Sys.Date(), "_MLNEfocal_quantiles", paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 1 (open pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(1, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(1, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	# num alleles per locus (only for adults and juvs)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = genowide$sizeclass == 4 | genowide$sizeclass == 1
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
	a = intersect(which(genowide$sizeclass == 4), sites)
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
	a = intersect(which(genowide$sizeclass == 1), sites)
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
}



###############################################
## Output to MLNE2: All samples (upper quartile vs. lower quartile, some pops as focal, certain others as source)

# sets of 1 w/ surrounding two as source
focalpops = list(b=8, c=9, d=10, e=11)

sourcepops = list(b=c(7,9), c=c(8,10), d=c(9,11), e=c(10,13))

prefix="quantilesSurrTwo"

NeMax = 10000

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75), na.rm=TRUE)

first = genowide$Size <= quants[1]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[fourth_inclusive] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]


# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
for(rep in 1:length(focalpops)){	
	focal = focalpops[[rep]]	
	sites = matchall(focal, genowide$SiteNum)
	a = intersect(which(genowide$sizeclass == 4), sites)
	b = intersect(which(genowide$sizeclass == 1), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}


# Write all groups to files
for(rep in 1:length(focalpops)){
	
	focal = focalpops[[rep]]
	source = sourcepops[[rep]]

	# Write to file
	mlnefile = file(paste("../Datafiles/MLNE/PRBI_", Sys.Date(), "_MLNEfocal_", prefix, paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 1 (open pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(1, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(1, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	sites = matchall(focal, genowide$SiteNum)
	source_sites = matchall(source, genowide$SiteNum)

	# num alleles per locus (only for adults and juvs)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = intersect(which(genowide$sizeclass==1 | genowide$sizeclass == 4), c(sites, source_sites))
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
	
	# For adults in focal pops:
	# Numbers of copies of each allele, at each locus, from each sample of the focal population
	cols = grep("A1.", names(genowide))
	a = intersect(which(genowide$sizeclass==4), sites)
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
	a = intersect(which(genowide$sizeclass==1), sites)
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

	# For source pop (adults and juvs)
	cols = grep("A1.", names(genowide))
	a = intersect(which(genowide$sizeclass==4 | genowide$sizeclass==1), source_sites)
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
}
















###############################
#### LEGACY CODE FROM APCL ####
###############################


########################
# Format for GenAlEx 6

	# Write to file
	genalexfile = file(paste("Aclarkii_GenAlEx_", Sys.Date(), ".csv", sep=""))
	open(genalexfile, "w")

	# Header: num loci, num samples, num pops, # indiv in each pop, # regions, # indivs in each region
	outline = c(14, dim(genowide)[1], length(unique(genowide$SiteNum)), table(genowide$SiteNum), length(unique(genowide$Region)), table(genowide$Region))
	cat(outline, file=genalexfile, sep=",", append=F)
	cat("\n", file=genalexfile, append=T)
	
	outline = c("Aclarkii 3 regions", "", "", names(table(genowide$SiteNum)), "", names(table(genowide$Region)))
	cat(outline, file=genalexfile, sep=",", append=T)
	cat("\n", file=genalexfile, append=T)

	i = grep("A1", names(genowide))
	j = gsub("A1.", "", names(genowide)[i], fixed=T)
	j = paste(j, collapse=",,")
	outline = paste("Sample no.,Pop,",j, ",,,Lat,Long",sep="") 
	cat(outline, file=genalexfile, append=T)
	cat("\n", file=genalexfile, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[,c(grep("Sample",widenames),grep("SiteNum",widenames),grep("A[[:digit:]]",widenames),grep("blank",widenames),grep("lat",widenames), grep("long",widenames))]
	write.table(outfile, file=genalexfile, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(genalexfile)	
	
	
	
###############################################
## Divide by size classes, output to GenAlEx

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)

dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)


hist(genowide$Size)
quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]

sum(first)
sum(second)
sum(third)
sum(fourth)

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4

# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

	# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
	widenames = names(genowide)
	i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
	genowide = genowide[,-i]
	
	
	# Write to file
	genalexfile = file(paste("Aclarkii_sizes_", Sys.Date(), "_GX.csv", sep=""))
	open(genalexfile, "w")

	# Header: num loci, num samples, num pops, # indiv in each pop, # regions, # indivs in each region
	outline = c(13, dim(genowide)[1], length(unique(genowide$sizeclass)), table(genowide$sizeclass))
	cat(outline, file=genalexfile, sep=",", append=F)
	cat("\n", file=genalexfile, append=T)
	
	outline = c("Aclarkii sizeclasses", "", "", names(table(genowide$sizeclass)))
	cat(outline, file=genalexfile, sep=",", append=T)
	cat("\n", file=genalexfile, append=T)

	i = grep("A1", names(genowide))
	j = gsub("A1.", "", names(genowide)[i], fixed=T)
	j = paste(j, collapse=",,")
	outline = paste("Sample no.,Pop,",j, ",,,Lat,Long",sep="") 
	cat(outline, file=genalexfile, append=T)
	cat("\n", file=genalexfile, append=T)

	# Data: Sample, Pop, cols of genotype data, blank col, Lat, Long
	genowide$blank = ""
	widenames = names(genowide)
	outfile = genowide[,c(grep("Sample",widenames),grep("sizeclass",widenames),grep("A[[:digit:]]",widenames),grep("blank",widenames),grep("lat",widenames), grep("long",widenames))]
	write.table(outfile, file=genalexfile, append=T, quote=F,row.names=F, col.names=F, sep=",", na="0")
	
	
	close(genalexfile)
	
	



###############################################
## Output to MLNE2: Focal pops as closed (TopTwo on anemone >= 8cm vs. lower quartile)

# sets of 4
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27))
prefix="TopTwo"

# sets of 4, 3, 2, and 1 (for testing only)
focalpops = list(a = c(7,8,9,10), b = c(8,9,10,11), c = c(9,10,11,13), d = c(10,11,13,14), e = c(11,13,14,15), f = c(13,14,15,16), g = c(14,15,16,17), h = c(18,19,20,22), i = c(19,20,22,23), j = c(20,22,23,24), k = c(22,23,24,25), l = c(23,24,25,27), m = c(7,8,9), n = c(8,9,10), o=c(9,10,11), p=c(10,11,13), q=c(11,13,14), r=c(13,14,15), s=c(14,15,16), t=c(15,16,17), u=c(18,19,20), v=c(19,20,22), w=c(20,22,23), x=c(22,23,24), y=c(23,24,25), z=c(24,25,27), aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27), qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)

# sets of 1,2, and 3, all using pop 25 (determined from sets above)
focalpops = list(a = 25, b = c(24,25), c=c(23,24,25))
prefix="TopTwo"

# sets of 2
focalpops = list(aa=c(7,8), bb=c(8,9), cc=c(9,10), dd=c(10,11), ee=c(11,13), ff=c(13,14), gg=c(14,15), hh=c(15,16), ii=c(16,17), jj=c(18,19), kk=c(19,20), ll=c(20,22), mm=c(22,23), nn=c(23,24), oo=c(24,25), pp=c(25,27))
prefix="TopTwo"

# sets of 1
focalpops = list(qq=7, rr=8, ss=9, tt=10, uu=11, vv=13, ww=14, xx=15, yy=16, zz=17, aaa=18, bbb=19, ccc=20, ddd=22, eee=23, fff=24, ggg=25, hhh=27)
prefix="TopTwo"


NeMax = 10000

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	
	
	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
for(rep in 1:length(focalpops)){	
	focal = focalpops[[rep]]	
	sites = matchall(focal, genowide$SiteNum)
	a = intersect(which(genowide$TopTwo & genowide$Size >= 8), sites)
	b = intersect(which(genowide$sizeclass == 1), sites)
	print(paste("Pops ", paste(focal, collapse=","), ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"

# Write all groups to files
for(rep in 1:length(focalpops)){
	
	focal = focalpops[[rep]]

	# Write to file
	mlnefile = file(paste("MLNE/Aclarkii_", Sys.Date(), "_MLNEclosed", prefix, paste(focal, collapse=""), ".csv", sep=""))
	open(mlnefile, "w")

	# Header: 0 (closed pop), then max Ne, then screen indicator, num cpus (0 for all), then num loci, 
	cat(0, file=mlnefile, sep=",", append=F)
	cat("\n", file= mlnefile, append=T)
	cat(NeMax, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)
	cat(2, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T) 
	cat(1, file=mlnefile, sep=",", append=T) # number of CPUs
	cat("\n", file= mlnefile, append=T)
	cat(13, file=mlnefile, sep=",", append=T)
	cat("\n", file= mlnefile, append=T)

	sites = matchall(focal, genowide$SiteNum)

	# num alleles per locus (only for adults and juvs in focal pop)
	outline = numeric()
	cols = grep("A1.", names(genowide))
	alleles = vector("list", 13)
	a = intersect(which(genowide$AdJuv == "Ad" | genowide$AdJuv == "Juv"), sites)
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
	
	close(mlnefile)	
}




###############################################
## Output for Geneland

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	

gnlndfile =(paste("Aclarkii_", Sys.Date(), "_Gnlndgen.csv", sep=""))
geofile =(paste("Aclarkii_", Sys.Date(), "_Gnlndgeo.csv", sep=""))

# Write out genetic data
write.table(subset(genowide, select=c(-Sample)), gnlndfile, sep=" ", row.names=F, col.names=F)
	
	# add lat/long info
	i = genowide$Sample
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(genowide, subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), all.x=T, by="Sample", sort=F)
	dim(genowide)
	
	identical(i, genowide$Sample) # must be TRUE: same order of samples
	
	
# Write out geographic data in Lambert coordinates
## load packages 
require(mapproj) 
require(maps) 
## check 
plot(genowide$long, genowide$lat,type="n",xlab="Lon",ylab="Lat",asp=1) 
points(genowide$long, genowide$lat,col=2) 
map(resol=0,add=TRUE) 
## convert (Lon,Lat) coordinates into Lambert 
mapproj.res <- mapproject(x=genowide$long, y=genowide$lat, projection="lambert", 
param=c(min(genowide$lat),max(genowide$lat))) 
## save planar coordinates as a two-column matrix 
coord.lamb <- cbind(mapproj.res$x,mapproj.res$y) 
colnames(coord.lamb) = c("X", "Y")
write.table(coord.lamb, geofile, row.names=F, col.names=F)



###############################################
## Output for structure

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	


	# Write to file
	structfile =file(paste("Aclarkii_", Sys.Date(), "_struct.csv", sep=""))
	open(structfile, "w")

	# Header: loci names
	i = grep("A1", names(genowide))
	j = gsub("A1.", "", names(genowide)[i], fixed=T)
	j = paste(j, collapse="  ")
	j = paste(" ",j,sep="")
	cat(j, file= structfile, sep=" ", append=F)
	cat("\n", file= structfile, append=T)

	# Data: Sample, cols of genotype data
	widenames = names(genowide)
	outfile = genowide[,c(grep("Sample",widenames),grep("A[[:digit:]]",widenames))]
	write.table(outfile, file= structfile, append=T, quote=F,row.names=F, col.names=F, sep=" ", na="-99")
	
	
	close(structfile)


###############################################
## Output for Alleles in Space (Miller 2005 J Hered)

geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	

aisfile =(paste("Aclarkii_", Sys.Date(), "_AISgen.csv", sep=""))
geofile =(paste("Aclarkii_", Sys.Date(), "_AISgeo.csv", sep=""))

# Write out genetic data
write.table(13,aisfile, sep="", row.names=F, col.names=F, append=F) # num loci

# collapse locus alleles with \
genowide[is.na(genowide)] = 0
cols = grep("A1.", names(genowide))
out = as.character(genowide$Sample)
out = cbind(out, paste(genowide[,cols[1]], genowide[,(cols[1]+1)], sep="\\"))
for(i in cols[2:length(cols)]){
	out = cbind(out, paste(genowide[,i], genowide[,(i+1)], sep="\\"))
}
write.table(out, aisfile, sep=", ", row.names=F, col.names=F, append=T, quote=F)
write.table(";", aisfile, sep=", ", row.names=F, col.names=F, append=T, quote=F)
	
	# add lat/long info
	i = genowide$Sample
	locations = read.csv("../../Surveys/Collections2008-11-21.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(genowide, subset(locations, select=c(Sample,SurveyNum,Size,lat,long)), all.x=T, by="Sample", sort=F)
	dim(genowide)
	
	identical(i, genowide$Sample) # must be TRUE: same order of samples
	
	
# Write out geographic data in Lambert coordinates
## load packages 
require(mapproj) 
require(maps) 
## check 
plot(genowide$long, genowide$lat,type="n",xlab="Lon",ylab="Lat",asp=1) 
points(genowide$long, genowide$lat,col=2) 
map(resol=0,add=TRUE) 
## convert (Lon,Lat) coordinates into Lambert 
mapproj.res <- mapproject(x=genowide$long, y=genowide$lat, projection="lambert", 
param=c(min(genowide$lat),max(genowide$lat))) 
## save planar coordinates as a two-column matrix 
offset=1 # add an arbitrary offset so that all coords are positive
coord.lamb <- cbind(as.character(genowide$Sample), (mapproj.res$x+offset),(mapproj.res$y+offset)) 
write.table(coord.lamb, geofile, row.names=F, col.names=F, sep=", ", quote=F, append=F)
write.table(";", geofile, sep=", ", row.names=F, col.names=F, append=T, quote=F)




###############################################
## TMVP: Output juvs and adults from one focal pop: TopTwo on anemone >= 8cm vs. lower quartile)
## Uses Genepop format with two populations. Run through Formatomatic to get to TMVP

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)
source("../Analysis/matchall.R")	
	
# reshape to wide
dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames
dim(genowide)

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]
	
	
	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)
	dim(genowide)

quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))

first = genowide$Size <= quants[1]
second = genowide$Size > quants[1] & genowide$Size <= quants[2]
third = genowide$Size > quants[2] & genowide$Size <= quants[3]
fourth = genowide$Size > quants[3]
fourth_inclusive = genowide$Size >= quants[3] # the largest fish

genowide$sizeclass = NA
genowide$sizeclass[first] = 1
genowide$sizeclass[second] = 2
genowide$sizeclass[third] = 3
genowide$sizeclass[fourth] = 4 # the largest fish


# order by sizeclass
i = order(genowide$sizeclass)
genowide = genowide[i,]

genowide$AdJuv = "" 
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad"
genowide$AdJuv[genowide$sizeclass == 1] = "Juv"


# How many breeding adults (TopTwo and > 8cm) and juvs (lowest quartile) in each population group?
groups =c("Ad", "Juv")
pops = sort(unique(genowide$SiteNum))
for(j in 1:length(pops)){	
	a = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Ad")
	b = which(genowide$SiteNum==pops[j] & genowide$AdJuv == "Juv")
	print(paste("Pops ", pops[j], ": Num Adults: ", length(a), "   Num Juvs: ", length(b), sep=""))
}

# Write out focal pop
pop = 25
gens = c(0,1) # generation times of the groups
cols = grep("A1.", names(genowide))
loci = unlist(strsplit(names(genowide)[cols], "A1.", fixed=T))[seq(2,26, by=2)]
	
	# Write to file
	name = paste("Aclarkii_", Sys.Date(), "_TVMP", pop, ".gen", sep="")
	file = file(name)
	open(file, "w")

	# Header
	cat(paste("Title line: Aclarkii pops, ", pop, " Ads and Juvs", sep=""), file= file, sep=",", append=F)
	cat("\n", file= file, append=T)
	cat(paste(loci, collapse="\n"), file= file, append=T)
	cat("\n", file= file, append=T) 
	for(i in 1:length(groups)){
		k = which(genowide$SiteNum==pop & genowide$AdJuv == groups[i])
		if(length(k)>0){
			cat("Pop", file= file, sep="", append=T)
			cat("\n", file= file, append=T)
		
			# collapse locus alleles with ""
			out = paste("     ", rep(gens[i], length(k)), ",", sep="") # name each indiv according to population name
			out = cbind(out, paste(formatC(genowide[k,cols[1]], width=3, flag="0"), formatC(genowide[k,(cols[1]+1)], width=3, flag="0"), sep=""))
			for(c in cols[2:length(cols)]){
				out = cbind(out, paste(formatC(genowide[k,c], width=3, flag="0"), formatC(genowide[k,(c+1)], width=3, flag="0"), sep=""))
			}
			out[out==" NA NA"] = "000000"
			for(c in 1:dim(out)[1]){
				cat(paste(out[c,], collapse=" "), file=file, append=T)
				cat("\n", file= file, append=T)
			}	
		}
	}
	close(file)	
}



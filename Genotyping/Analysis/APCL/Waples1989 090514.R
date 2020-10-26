setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Genotypes")
source("../Analysis/matchall.R")	

###########################
# Read in and format data
geno = read.csv("Aclarkii_genotypes_2009-03-13.csv", row.names=1)

dropnames = names(geno)
dropnames = c(dropnames[grep("A[[:digit:]][[:punct:]]", dropnames)], "error", "numgenos")
genowide = reshape(geno, direction="wide", v.names = c("A1consens", "A2consens"), timevar = "Marker", idvar = "Sample", drop=dropnames)
widenames = names(genowide)
widenames = gsub("consens", "", widenames)
names(genowide) = widenames

# remove ACH_A7, ACH_B9 and APR_Cf42 because out of HWE, high error or low completeness (3/15/09)
widenames = names(genowide)
i = c(grep("ACH_A7", widenames), grep("ACH_B9", widenames), grep("APR_Cf42", widenames))
genowide = genowide[,-i]

# Drop APCL285, 286 since identical
i = genowide$Sample == "APCL285" | genowide$Sample == "APCL286"
genowide = genowide[-i,]
	
	
	# add population and lat/long info
	surveys = read.csv("../../Surveys/GPSSurveys2009-01-08.csv")
	locations = read.csv("../../Surveys/Collections2009-03-26.csv")
	locations$Sample = paste(locations$Spp, locations$ID, sep="")
	genowide = merge(subset(locations, select=c(Sample,SurveyNum,Size,lat,long, TopTwo)), genowide, all.y=T, by="Sample")
	genowide = merge(subset(surveys, select=c(SurveyNum,SiteNum,Name,Region,Municipality)), genowide, all.y=T)

genowide$AdJuv = "" 
genowide$sizeclass = 0
quants = quantile(genowide$Size, probs = c(0.25, 0.5, 0.75))
genowide$sizeclass[genowide$Size <= quants[1]] = 1
genowide$sizeclass[genowide$Size > quants[3]] = 4
genowide$AdJuv[genowide$TopTwo & genowide$Size >= 8] = "Ad" # an adult according to largest pair and >8cm
genowide$AdJuv[genowide$Size <= quants[1]] = "Juv"
i = order(genowide$sizeclass)
genowide = genowide[i,]


#############################
## Calculate Ne

# test the waples89 function against my Excel spreadsheet calcs: 1000 (737, 1302)
genojuv = genowide[genowide$sizeclass == 1, ]
genoad =  genowide[genowide$sizeclass == 4, ]

waples89(genojuv, genoad,1)

# test the function with pooling alleles < 0.02
waples89(genojuv, genoad, 1, pool=0.02)

# try calc using Ad/Juv distinction
genojuv = genowide[genowide$AdJuv == "Juv", ]
genoad =  genowide[genowide$AdJuv == "Ad", ]

waples89(genojuv, genoad,1)
waples89(genojuv, genoad,1/3) # try 1/3 generation between samples


# Ad/Juv and pooling < 0.02
genojuv = genowide[genowide$AdJuv == "Juv", ]
genoad =  genowide[genowide$AdJuv == "Ad", ]
waples89(genojuv, genoad, 1, pool=0.02)


# single site: doesn't work
sites = sort(unique(genowide$SiteNum))
len = length(sites)
onesite = data.frame(site = sites, ne = numeric(len), L95 = numeric(len), U95 = numeric(len))

for(i in 1:len){
	genojuv = genowide[genowide$AdJuv=="Juv" & genowide$SiteNum==sites[i],]
	genoad = genowide[genowide$AdJuv=="Ad" & genowide$SiteNum==sites[i],]
	print(paste(dim(genojuv)[1],",",dim(genoad)[1]))
	out = waples89(genojuv, genoad,2)
	#print(out)
	onesite[i,2:4] = out
}
round(onesite)


# four sites: doesn't work well
for(i in sort(unique(genowide$SiteNum))){
	genojuv = genowide[genowide$AdJuv=="Juv" & (genowide$SiteNum==i | genowide$SiteNum==(i+1)| genowide$SiteNum==(i+2)| genowide$SiteNum==(i+3)),]
	genoad = genowide[genowide$AdJuv=="Ad" & (genowide$SiteNum==i | genowide$SiteNum==(i+1)| genowide$SiteNum==(i+2)| genowide$SiteNum==(i+3)),]
	out = c(paste("Site", i), round(waples89(genojuv, genoad,1)))
	print(out)
}

# eight sites: works a little better
for(i in sort(unique(genowide$SiteNum))){
	genojuv = genowide[genowide$AdJuv=="Juv" & (genowide$SiteNum==i | genowide$SiteNum==(i+1)| genowide$SiteNum==(i+2)| genowide$SiteNum==(i+3)| genowide$SiteNum==(i+4)| genowide$SiteNum==(i+5)| genowide$SiteNum==(i+6)| genowide$SiteNum==(i+7)),]
	genoad = genowide[genowide$AdJuv=="Ad" & (genowide$SiteNum==i | genowide$SiteNum==(i+1)| genowide$SiteNum==(i+2)| genowide$SiteNum==(i+3)| genowide$SiteNum==(i+4)| genowide$SiteNum==(i+5)| genowide$SiteNum==(i+6)| genowide$SiteNum==(i+7)),]
	out = c(paste("Site", i), round(waples89(genojuv, genoad,1)))
	print(out)
}

# Cebu vs. Leyte: only works on Leyte (larger sample size....)
genojuv = genowide[genowide$AdJuv=="Juv" & genowide$Region=="Cebu",]
genoad = genowide[genowide$AdJuv=="Ad" & genowide$Region=="Cebu",]
dim(genojuv)
dim(genoad)
waples89(genojuv, genoad,1)
waples89(genojuv, genoad,1, pool=0.05)

genojuv = genowide[genowide$AdJuv=="Juv" & genowide$Region=="Leyte",]
genoad = genowide[genowide$AdJuv=="Ad" & genowide$Region=="Leyte",]
dim(genojuv)
dim(genoad)
waples89(genojuv, genoad,1)
waples89(genojuv, genoad,1, pool=0.05)


#####################################################

waples89 <- function(genojuv, genoad, g, pool = NA){  # give it a subset of the genowide table corresponding to juvs and adults
	 # g = number of generations between samples
	 # pool = threshhold below which to pool alleles

	#################################
	## Set up allele frequency table (af) and sample size vector (S)
	# rows for each allele of each locus, cols for frequency in juvs and adults

	cols = grep("A1.", names(genowide))
	af = data.frame(locus = character(0), allele = numeric(0), ad = numeric(0), juv = numeric(0))
	S = data.frame(allele = character(0), juv = numeric(0), ad = numeric(0))
	for(i in 1:length(cols)){ # for each locus
		# juvs
		jals = c(genojuv[,cols[i]], genojuv[,cols[i]+1])
		jals = jals[!is.na(jals)]
		freqs = unlist(lapply(split(jals,f=jals),length))/length(jals)
		juv = data.frame(allele = as.numeric(names(freqs)), juv = freqs)
	
		# adults
		aals = c(genoad[,cols[i]], genoad[,cols[i]+1])
		aals = aals[!is.na(aals)]
		freqs = unlist(lapply(split(aals,f=aals),length))/length(aals)
		ad = data.frame(allele = as.numeric(names(freqs)), ad = freqs)

		temp = merge(juv, ad, all.x=T, all.y=T, sort=T)
		temp$locus = strsplit(names(genowide)[cols[i]], ".", fixed=T)[[1]][2]
	
		af = rbind(af, temp)
		S = rbind(S, data.frame(locus = temp$locus[1], juv = length(jals)/2, ad = length(aals)/2))
	}
	af[is.na(af)] = 0
	dim(af)
	
	# pool alleles if required
	if(!is.na(pool)){
		loc = unique(af$locus)
		for(i in 1:length(loc)){
			k = which(af$juv < pool & af$ad < pool & af$locus == loc[i])
			if(length(k)>0){
				temp = data.frame(locus = loc[i], allele = "pool", juv = sum(af$juv[k]), ad = sum(af$ad[k]))
				af = rbind(af[-k,], temp)
			}
		}
	}
	dim(af)

	#### Calculations
	af$xiyi = (af$ad-af$juv)^2/((af$ad+af$juv)/2-af$ad*af$juv)
	Kj = unlist(lapply(split(af$locus, f=af$locus), length))
	F = unlist(lapply(split(af$xiyi, f=af$locus), sum))/Kj
	n = sum(Kj-1) # number of independent alleles
	Fbar = sum(Kj*F)/sum(Kj)
	Fbar = c(Fbar, n*Fbar/qchisq(0.025, df=n), n*Fbar/qchisq(0.975, df=n)) # upper and lower bounds on Fbar
	
	# correction factor 
	Sjuv = sum(Kj-1)/sum((Kj-1)/S$juv)
	Sad = sum(Kj-1)/sum((Kj-1)/S$ad)
	corr = 1/(2*Sjuv) + 1/(2*Sad)

	Fprime = Fbar-corr
	invN = 0 # 1/N, can be set to about 0
	
	Ne = g/(2*(Fprime+invN))
	names(Ne) = c("mu", "L95", "U95")
	return(Ne)
}



#######################################
## Krimbas & Tsakas 1971 example data

af = data.frame(locus = c(rep("A", 17), rep("B", 13)), allele=c("00", "0", "01", "1", "112", "12", "2", "23", "3", "34", "4", "5", "6", "7", "78", "8", "s", "0", "1", "2", "3", "4", "45", "5", "56", "566", "6", "7", "8", "s"), juv=c(.00105, .01266, .00105, .13947, 0, .04865, .43558, .00847, .05533, 0, .02456, .02456, .00844, .08713, 0, .03434, .11871, 0, .00107, .04807, 0.02046, .62203, .03198, .03693, .01279, .00426, .07533, .09987, .00320, .04410), ad=c(0, .00801, 0, .11390, .0016, .04079, .45404, .01616, .03424, .00160, .03590, .01616, .02434, .06802, .00321, .01603, .16600, .00177, .00712, .04739, .01615, .63713, .01068, 0.05677, .01068, 0, 0.06625, .10319, 0, .04287))

S = data.frame(locus = c("A", "B"), juv=c(474, 469), ad = c(312, 281))
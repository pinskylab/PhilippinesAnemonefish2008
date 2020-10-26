setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys/Genotype Analysis/FaMoz")

#famoz = readLines("Famoz 090416.txt", n=-1)
famoz = readLines("090421/Famoz 090421.txt", n=-1)
ids = read.table("Aclarkii_create_2009-04-15-FAMOZ-ID TABLE.txt", sep="\t", header=TRUE)


#famoz = readLines("090423/Famoz 090423.txt", n=-1)
#ids=read.table("Aclarkii_create_2009-04-23no18or7-FAMOZ-ID TABLE.txt", sep="\t", header=TRUE)

#famoz = readLines("090424/Famoz 090424.txt", n=-1)
#ids=read.table("Aclarkii_create_2009-04-24_pops2223-FAMOZ-ID TABLE.txt", sep="\t", header=TRUE)


loc = read.csv("Collections2009-03-26.csv")
loc$name = paste(loc$Spp, loc$ID, sep="")
names(loc)[grep("lat", names(loc))] = "latnew"
names(loc)[grep("long", names(loc))] = "longnew"

#090416 threshholds, err=0.001
#threshone = 3.91 # threshold LOD score for assigning a parent
#threshpair = 10.72 # threshold LOD score for assigning a parent pair

#extreme thresholds
#threshone = 7 # threshold LOD score for assigning a parent
#threshpair = 15 # threshold LOD score for assigning a parent pair

# 090420 thresholds, err=0
threshone = 3.38 # threshold LOD score for assigning a parent
threshpair = 9 # threshold LOD score for assigning a parent pair


# annotate ids as parent or juv
len = length(ids$FAMOZ.ID)
i = ids$FAMOZ.ID[2:len] - ids$FAMOZ.ID[1:(len-1)] # find where id drops back to 1
i = which(i<0)
ids$adjuv = 0  # adults
ids$adjuv[(i+1):len] = 1 # juveniles

kids = ids[ids$adjuv == 1, c(1,2)]
kids$numpars = 0 # number of likely parents
kids$numpairs = 0 # number of likely parent pairs
kids$p1 = "" # parent one name
kids$p2 = ""
kids$p3 = ""
kids$p4 = ""
kids$p1lod = 0 # parent one LOD score
kids$p2lod = 0
kids$p3lod = 0
kids$p4lod = 0
kids$pp11 = "" # name of parent one in parent pair one
kids$pp12 = "" # name of parent two in parent pair one
kids$pp21 = ""
kids$pp22 = ""
kids$pp31 = ""
kids$pp32 = ""
kids$pp41 = ""
kids$pp42 = ""
kids$pp1lod = 0 # parent pair 1 LOD score
kids$pp2lod = 0
kids$pp3lod = 0
kids$pp4lod = 0


# read in parents
i = grep("***Most likely parents and parent pair for each offspring ***", famoz, fixed=TRUE)
if(length(i)>1){ print("To many parentage assignments in this file!") }

kidlines = grep("kid", famoz, fixed=TRUE)
length(kidlines) # should equal the number of kids

for(i in 1:length(kidlines)){
	j = kidlines[i]
	thiskid = as.numeric(unlist(strsplit(famoz[j], " ", fixed=TRUE))[3])
	if(i != thiskid) {print(paste("Why isn't thiskid == i?", i))}
	
	outline = which(kids$FAMOZ.ID == thiskid) # get the line index into kids to write to

	j = j+1 # skip to num of parents
	numpars = as.numeric(unlist(strsplit(famoz[j], " ", fixed=TRUE))[2])
	kids$numpars[outline] = numpars
	j = j+1
	if(numpars > 0){ # if there are parents to read
		for(k in 1:min(numpars,4)){ # read in each likely parent, up to 4
			temp = unlist(strsplit(famoz[j], "\t"))
			namecol = paste("p", k, sep="")
			lodcol = paste("p", k, "lod", sep="")
			par = as.numeric(temp[2]) # name of parent
			par = as.character(ids$RAW.ID[which(ids$adjuv == 0 & ids$FAMOZ.ID == par)])
			kids[[namecol]][outline] = par
			kids[[lodcol]][outline] = as.numeric(temp[3]) # lod score
			j = j+1 # move to next parent or to parent pairs
		}
	}
	numpairs = as.numeric(unlist(strsplit(famoz[j], " ", fixed=TRUE))[2])
	kids$numpairs[outline] = numpairs
	j = j+1
	if(numpairs > 0){ # if there are parent pairs to read
		for(k in 1:min(numpairs,4)){ # read in each likely pair, up to 4
			temp = unlist(strsplit(famoz[j], "\t"))
			namecol1 = paste("pp", k, "1", sep="")
			namecol2 = paste("pp", k, "2", sep="")
			lodcol = paste("pp", k, "lod", sep="")
			pars = unlist(strsplit(gsub("[[:punct:]]", "", temp[2]), " ", fixed=T))[2:3]
			par = as.character(ids$RAW.ID[which(ids$adjuv == 0 & ids$FAMOZ.ID == as.numeric(pars[1]))])
			kids[[namecol1]][outline] = par # name of parent 1
			par = as.character(ids$RAW.ID[which(ids$adjuv == 0 & ids$FAMOZ.ID == as.numeric(pars[2]))])
			kids[[namecol2]][outline] = par # name of parent 2
			kids[[lodcol]][outline] = as.numeric(temp[3]) # lod score
			j = j+1 # move to next parent or to parent pairs
		}
	}
}

# take a subset of the most likely
kids2 = subset(kids, select=c("RAW.ID", "FAMOZ.ID", "p1", "p1lod", "pp11", "pp12", "pp1lod"))


# erase parents and pairs with lod < threshold 
i = kids2$p1lod < threshone
	kids2$p1[i] = ""
	kids2$p1lod[i] = 0
	kids2$pp1lod[i] = 0
i = kids2$pp1lod < threshpair
	kids2$pp11[i] = ""
	kids2$pp12[i] = ""
	kids2$pp1lod[i] = 0

# delete the single parents if there's a pair
i = kids2$pp1lod >= threshpair
	kids2$p1[i] = ""
	kids2$p1lod[i] = ""


# How many single parents and pairs?
print(paste("Number of kids:", length(kids2$RAW.ID)))
i = sum(kids2$p1lod>0)
print(paste("Number of single parents:", i))
i = sum(kids2$pp1lod>0)
print(paste("Number of parents pairs:", i))


# add lat and long
kids2 = merge(kids2, subset(loc, select=c(name, latnew, longnew)), by.x="RAW.ID", by.y="name", all.x=TRUE)
names(kids2)[grep("latnew", names(kids2))] = "lat.kid"
names(kids2)[grep("longnew", names(kids2))] = "long.kid"

kids2 = merge(kids2, subset(loc, select=c(name, latnew, longnew)), by.x="p1", by.y="name", all.x=TRUE)
names(kids2)[grep("latnew", names(kids2))] = "lat.p1"
names(kids2)[grep("longnew", names(kids2))] = "long.p1"

kids2 = merge(kids2, subset(loc, select=c(name, latnew, longnew)), by.x="pp11", by.y="name", all.x=TRUE)
names(kids2)[grep("latnew", names(kids2))] = "lat.pp11"
names(kids2)[grep("longnew", names(kids2))] = "long.pp11"

kids2 = merge(kids2, subset(loc, select=c(name, latnew, longnew)), by.x="pp12", by.y="name", all.x=TRUE)
names(kids2)[grep("latnew", names(kids2))] = "lat.pp12"
names(kids2)[grep("longnew", names(kids2))] = "long.pp12"

# For parent pairs: are parents on the same anemone?
k = which(kids2$pp1lod > 0)
i = sum(kids2$lat.pp11[k] == kids2$lat.pp12[k])
j = sum(kids2$long.pp11[k] == kids2$long.pp12[k])
print(paste("Parents at same lat: ", round(i/length(k)*100),"%", sep=""))
print(paste("Parents at same long: ", round(j/length(k)*100),"%", sep=""))

# Calculate distance
require(argosfilter)
kids2$dist1km = NA # distance to parent 1
kids2$dist2km = NA # distance to parent 2
kids2$distppkm = NA # distance between parents
k = which(kids2$p1lod>0)
for(i in k){
	kids2$dist1km[i]=distance(lat1=kids2$lat.kid[i], lat2=kids2$lat.p1[i], lon1=kids2$long.kid[i], lon2=kids2$long.p1[i])
}
k = which(kids2$pp1lod>0)
for(i in k){
	kids2$dist1km[i]=distance(lat1=kids2$lat.kid[i], lat2=kids2$lat.pp11[i], lon1=kids2$long.kid[i], lon2=kids2$long.pp11[i])
	kids2$dist2km[i]=distance(lat1=kids2$lat.kid[i], lat2=kids2$lat.pp12[i], lon1=kids2$long.kid[i], lon2=kids2$long.pp12[i])
	kids2$distppkm[i]=distance(lat1=kids2$lat.pp11[i], lat2=kids2$lat.pp12[i], lon1=kids2$long.pp11[i], lon2=kids2$long.pp12[i])
}



## Graph parent-offspring and parent-parent distances
par(mfrow=c(1,2))
hist(c(kids2$dist1km,kids2$dist2km), breaks=20, col="grey30", xlab="Distance between kid and parents (km)", xlim=c(0,250), main="FaMoz")

hist(kids2$distppkm, breaks=20, col="grey40", xlab="Distance between parents (km)", xlim=c(0,200), main="FaMoz")




######### Check relatedness of parents and offspring
# Calculated as proportion of loci (NOT alleles) where kids and parents match at least one allele
geno = read.csv("Aclarkii_create_2009-04-15.csv")

pairs = which(kids2$p1lod>0 | kids2$pp1lod>0)
relatedness = numeric(0)
relatedness2 = numeric(0)
for(i in 1:length(pairs)){
	j = which(as.character(geno$Indiv)== kids2$RAW.ID[pairs[i]])
	if(length(j)>1){ print(paste("Too many kids found", i))}
	thiskid = geno[j,4:29]
	
	if(kids2$p1lod[pairs[i]]>0){ # for single parents
		j = which(as.character(geno$Indiv)== kids2$p1[pairs[i]])
		if(length(j)>1){ print(paste("Too many single parents found", i))}
		thisparent1 = geno[j,4:29]
	}
	
	if(kids2$pp1lod[pairs[i]]>0){ # for couples
		j = which(as.character(geno$Indiv)== kids2$pp11[pairs[i]])
		if(length(j)>1){ print(paste("Too many parent1 in pair found", i))}
		thisparent1 = geno[j,4:29]
		j = which(as.character(geno$Indiv)== kids2$pp12[pairs[i]])
		if(length(j)>1){ print(paste("Too many parent2 in pair found", i))}
		thisparent2 = geno[j,4:29]
	}

	# test parent1 against kid
	matches = 0
	for(j in seq(1,25,2)){
		kidg = thiskid[j:(j+1)]
		parg = thisparent1[j:(j+1)]
		a1 = match(kidg[1], parg) # only returns the first match
		a2 = match(kidg[2], parg) 
		# To calc proportion of loci that match
		r = max(sum(!is.na(a1)), sum(!is.na(a2)))		

		# Set a match if kid or parent is missing a genotype
		if((kidg[1]==0 & kidg[2]==0) | (parg[1]==0 & parg[2]==0)){
			r = 1
		} 

		# To calc proportion of alleles that match
		#r = sum(!is.na(a1)) + sum(!is.na(a2)) # otherwise, sum the number of matches (0, 1, or 2)
		# make sure both kid alleles didn't match to same parent allele (homozygote kid, hetero parent)
		#if((kidg[1] == kidg[2]) & (parg[1] != parg[2])){ 
		#	r = min(r,1)
		#}
		if(r<1){
			print(paste("No match for ", kids2$RAW.ID[pairs[i]], " at ", names(kidg)[1], sep=""))
		}		
		matches = matches + r
	}
	relatedness = c(relatedness, matches/26)

	# if there's a second parent (a couple)
	if(kids2$pp1lod[pairs[i]]>0){
		matches = 0
		for(j in seq(1,25,2)){
			kidg = thiskid[j:(j+1)]
			parg = thisparent2[j:(j+1)]
			a1 = match(kidg[1], parg) # only returns the first match
			a2 = match(kidg[2], parg) 
			# To calc proportion of loci that match
			r = max(sum(!is.na(a1)), sum(!is.na(a2)))		

			# Set a match if kid or parent is missing a genotype
			if((kidg[1]==0 & kidg[2]==0) | (parg[1]==0 & parg[2]==0)){
				r = 1
			}
	
			# To calc proportion of alleles that match
			#r = sum(!is.na(a1)) + sum(!is.na(a2)) # otherwise, sum the number of matches (0, 1, or 2)
			# make sure both kid alleles didn't match to same parent allele (homozygote kid, hetero parent)	
			#if((kidg[1] == kidg[2]) & (parg[1] != parg[2])){ 
			#	r = min(r,1)
			#}
			if(r<1){
				print(paste("No match for ", kids2$RAW.ID[pairs[i]], " at ", names(kidg)[1], sep=""))
			}		
			matches = matches + r
		}
		relatedness2 = c(relatedness2, matches/26)
	}	
}

r = c(relatedness, relatedness2)
summary(r)
hist(r)

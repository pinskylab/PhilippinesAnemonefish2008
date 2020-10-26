# Calculate relatedness for all possible pairs of fish
setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys/Genotype Analysis/FaMoz")

geno = read.csv("Aclarkii_create_2009-04-24.csv")

len = dim(geno)[1]
kd = character(0)
pr = character(0)
rl = numeric(0)
for(i in 1:(len-1)){
	cat(paste(i,"\n"))
	for(j in (i+1):len){
			thisparent1 = geno[i,4:29]
			thiskid = geno[j,4:29]
	
			# test parent1 against kid
			matches = 0
			locuscols = seq(1,25,2)
			for(k in locuscols){ # for each locus
				kidg = thiskid[k:(k+1)]
				parg = thisparent1[k:(k+1)]
				a1 = match(kidg[1], parg) # only returns the first match
				a2 = match(kidg[2], parg) 
				# To calc proportion of loci that match	
				r = max(sum(!is.na(a1)), sum(!is.na(a2)))	# did either allele match the parent?
	
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
		
				matches = matches + r
			}
			kd = c(kd, as.character(geno$Indiv[i]))
			pr = c(pr, as.character(geno$Indiv[j]))
			rl = c(rl, matches/length(locuscols))
			
	}
}
rel = data.frame(kid = kd, par = pr, rel = rl)

# Add cohort information
rel2 = merge(rel, subset(geno, select=c(Indiv, Cohort)), by.x = "kid", by.y="Indiv")
n = names(rel2)
i = grep("Cohort", n)
names(rel2)[i] = "kidcohort"
rel2 = merge(rel2, subset(geno, select=c(Indiv, Cohort)), by.x = "par", by.y="Indiv")
n = names(rel2)
i = grep("Cohort", n)
names(rel2)[i] = "parcohort"


# How well-related are all possible pairs of fish?
summary(rel$rel)
length(rel$rel)
sum(rel$rel==1)
sum(rel$rel==1)/length(rel$rel)
bk = (0:13+0.5)/13
hist(breaks=bk, rel$rel, xlab="% of loci with one or more shared alleles", main="Relatedness among Samples", col="grey")


# How well-related are all possible pairs of kids and parents?
i = which(rel2$kidcohort != rel2$parcohort) # a potential kid and parent were paired

summary(rel2$rel[i])
length(rel2$rel[i])
sum(rel2$rel[i]==1)
sum(rel2$rel[i]==1)/length(rel2$rel[i])
bk = (0:13+0.5)/13
hist(breaks=bk, rel2$rel[i], xlab="% of loci with one or more shared alleles", main="Relatedness among Kids and Parents", col="grey")

# How many unique kids and unique parents get matched?
i = which(rel2$kidcohort != rel2$parcohort)
j = rel2$rel[i]==1
sum(j)
length(k<-unique(geno$Indiv)) #376 indivs total
length(l<-unique(c(as.character(rel2$kid),as.character(rel2$par)))) # 376 indivs total
length(unique(c(as.character(rel2$kid[i]),as.character(rel2$par[i])))) # 376 indivs total
# 215 parents
length(k <- unique(c(as.character(rel2$kid[rel2$kidcohort==0]),
	as.character(rel2$par[rel2$parcohort==0]))))
# 161 kids
length(l <- unique(c(as.character(rel2$kid[rel2$kidcohort==1]),
	as.character(rel2$par[rel2$parcohort==1]))))
# 121 unique parents match a kid
length(unique(k<-c(as.character(rel2$kid[i][j][rel2$kidcohort[i][j]==0]),
	as.character(rel2$par[i][j][rel2$parcohort[i][j]==0]))))
	m = aggregate(k, by=list(k), length)
	summary(m[,2])
	n = hist(m[,2], breaks=seq(min(m[,2])-0.5, max(m[,2])+0.5, by=1), plot=FALSE)
	rbind(n$mids, n$counts)

# 102 unique kids match a parent
length(unique(k<-c(as.character(rel2$kid[i][j][rel2$kidcohort[i][j]==1]),
	as.character(rel2$par[i][j][rel2$parcohort[i][j]==1]))))
	m = aggregate(k, by=list(k), length)
	summary(m[,2])
	n = hist(m[,2], breaks=seq(min(m[,2])-0.5, max(m[,2])+0.5, by=1), plot=FALSE)
	rbind(n$mids, n$counts)

write.csv(rel2, paste("relatedness", Sys.Date(), ".csv", sep=""))

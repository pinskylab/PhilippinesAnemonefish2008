## Read Arlequin pairdist.sum file from jackknife over loci of CebuLeyte
## Write out .csvs of pairwise fsts

## NOT NEEDED: Use readArlJackknifes 100112.R

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Genotyping/Analysis/APCL/090507/Arlequin")
infile = readLines("jackknife/pairdist.sum") # read file

fileinds = grep(".arp", infile) # marks beginning of each Fst matrix
pvalinds = grep("p-Values", infile) # marks end of each Fst matrix
numfiles = length(fileinds)


fsts = vector("list", 0)


# read in fsts
for(i in 1:numfiles){
	rawfsts = strsplit(infile[(fileinds[i]+1):(pvalinds[i]-1)], "\t") # return fsts as a list
	numpops = length(rawfsts)
	thesefsts = matrix(data=NA, nrow=numpops, ncol=numpops)
	thesefsts[1,1] = 0
	for(j in 2:numpops){
		thesefsts[j,1:j] = as.numeric(rawfsts[[j-1]])
	}
	fsts[[i]] = thesefsts
	# add names
	splitname = unlist(strsplit(infile[fileinds[i]], ""))
	startname = max(grep("/", splitname))+1
	endname = grep("\"", splitname)
	endname = endname[min(which(endname>startname))]-5 # lop of .arp"
	names(fsts)[i] = paste(splitname[startname:endname], collapse="")

}




# Write results
for(i in 1:length(fsts)){
	outfile = paste("jackknife/Fst_summaries/", names(fsts)[i], "_fst.csv", sep="")
	write.csv(fsts[[i]], outfile)
}
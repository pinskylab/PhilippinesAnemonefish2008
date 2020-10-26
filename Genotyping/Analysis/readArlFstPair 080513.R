# Function to read in Arl Fsts and p-values from a directory of .res folders
# Reads Fsts and p-values from the pairwise FSTs between populations p1 and p2

readArlFstPair <- function(dirname, newdirname, p1, p2){
	dirs = list.files(dirname, pattern="^Set[[:digit:]]+$") # starts with "Set", has some numbers, and then ends (no .par, .sum, etc.)

	fsts = numeric()
	pvals = numeric()
	names = character()
	sets = character()

	# get our population numbers in order
	if(p1>p2){ 
		bigp = p1
		littlep = p2
	} else if(p2>p1){
		bigp = p2
		littlep = p1
	} else {
		warning("p1 and p2 must be different!")
	}

	for(i in 1:length(dirs)){ # iterate through all the "Set" directories
		cat(paste(dirs[i],"\n"))

		if(newdirname ==""){
			setdir = paste(dirname, "/", dirs[i], sep="")
		} else {
			setdir = paste(dirname, "/", dirs[i], "/", newdirname, sep="")
		}

		resultdirs = list.files(setdir, pattern=".res")
		
		for(j in 1:length(resultdirs)){ # iterate through all the .res folders in the Set directory
			temp = strsplit(resultdirs[j], split=".", fixed=T)[[1]][1]
			filename = paste(setdir, "/", resultdirs[j], "/", temp, ".htm", sep="")

			lines = readLines(filename) # read in this simulation's results
			fstline = grep("Population pairwise FSTs", lines)
			pvalline = grep("FST P values", lines)
			npops = pvalline - fstline - 9 # figure out how many populations this file has, from the size of the dist matrix

			if(bigp > npops | bigp< 1 | littlep < 1){
				warning("You picked a population that doesn't exist!")
			}

			if(length(fstline) != 0 & length(pvalline) != 0){
				line1 = fstline + 5 + bigp # get the correct line in the FST distance matrix
				line2 = pvalline + 5 + bigp # get the correct line in the FST distance matrix

				str1 = unlist(strsplit(lines[line1], "  "))
				str2 = unlist(strsplit(lines[line2], "  "))

				fsts = c(fsts, as.numeric(str1[6+littlep]))
				pvals = c(pvals, as.numeric(strsplit(str2[6+littlep], split="+", fixed=T)[[1]][1]))
				names = c(names, temp)
				sets = c(sets, dirs[i])

				
			} else {
					warning("Not an Arlequin .res file with Pairwise FSTs!")
			}
		}
	}

	out = data.frame(set = sets, name = names, fst = fsts, p = pvals)
	return(out)
	
}
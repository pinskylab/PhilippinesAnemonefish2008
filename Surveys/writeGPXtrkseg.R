# Write a trkseg to a GPX file

writeGPXtrkseg = function(outfilename, name, lat, long, time, startend){
	if(startend != "start" & startend != "end" & startend != "mid"){
		warning("Must specify start, mid, or end of track")
	}

	if(startend == "start"){
		cat("<trk>", file=outfilename, sep="\n", append=T)
		out = paste("<name>", name, "</name>", sep="")
		cat(out, file=outfilename, sep="\n", append=T)
	}


	cat("<trkseg>", file=outfilename, sep="\n", append=T)

	for(i in 1:length(lat)){
		out = paste("<trkpt lat=\"",lat[i], "\" lon=\"",long[i], "\">", sep="")
		cat(out, file=outfilename, sep="\n", append=T)
		cat("  <ele>0</ele>", file=outfilename, sep="\n", append=T)
		out = paste("  <time>",as.character(time[i]), "</time>", sep="")
		cat(out, file=outfilename, sep="\n", append=T)
		cat("</trkpt>", file=outfilename, sep="\n", append=T)
	}

	cat("</trkseg>", file=outfilename, sep="\n", append=T)
	
	if(startend == "end"){
		cat("</trk>", file=outfilename, sep="\n", append=T)
	}
}
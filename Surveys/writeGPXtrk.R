# Write a complete, continuous trk to a GPX file

writeGPXtrk = function(outfilename, name, lat, long, time){
	cat("<trk>", file=outfilename, sep="\n", append=T)
	out = paste("<name>", name, "</name>", sep="")
	cat(out, file=outfilename, sep="\n", append=T)
	cat("<trkseg>", file=outfilename, sep="\n", append=T)

	for(i in 1:length(lat)){
		out = paste("<trkpt lat=\"",lat[i], "\" lon=\"",long[i], "\">", sep="")
		cat(out, file=outfilename, sep="\n", append=T)
		cat("  <ele>0</ele>", file=outfilename, sep="\n", append=T)
		out = paste("  <time>",as.character(time[i]), "</time>", sep="")
		cat(out, file=outfilename, sep="\n", append=T)
		cat("</trkpt>", file=outfilename, sep="\n", append=T)
	}

	cat("</trkseg></trk>", file=outfilename, sep="\n", append=T)
	
}
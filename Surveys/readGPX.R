# Read a .gpx file and return a dataframe of time, lat, and long
# Reads ONLY from the "ACTIVE LOG"

readGPX = function(filename){
	infile = readLines(filename)
	lines = grep("<trkpt lat", infile, fixed=T)
	activestart = grep("<name>ACTIVE LOG</name>", infile, fixed=T)
	trksegend = grep("</trkseg>", infile, fixed=T)
	trksegend = trksegend[trksegend>activestart][1]
	lines = lines[lines>activestart & lines<trksegend] # Start with the beginning of the ACTIVE LOG and end with the end of its trkseg

	len = length(lines)
	
	time = character(len)
	year = numeric(len)
	month = numeric(len)
	day = numeric(len)
	hour = numeric(len)
	min = numeric(len)
	sec = numeric(len)
	lat = numeric(len)
	long = numeric(len)
	for(i in 1:len){
		thisline = lines[i]
		time[i] = strsplit(strsplit(infile[thisline+2],"<", fixed=T)[[1]][2], ">", fixed=T)[[1]][2]
		temp = strsplit(infile[thisline],"\"", fixed=T)[[1]]
		lat[i] = as.numeric(temp[2])
		long[i] = as.numeric(temp[4])

		temp = strsplit(time[i], "-", fixed=T)[[1]]
		year[i] = as.numeric(temp[1])
		month[i] = as.numeric(temp[2])
		temp = strsplit(temp[3],"T", fixed=T)[[1]]
		day[i] = as.numeric(temp[1])
		temp = strsplit(temp[2], ":", fixed=T)[[1]]
		hour[i] = as.numeric(temp[1])
		min[i] = as.numeric(temp[2])
		sec[i] = as.numeric(strsplit(temp[3],"Z",fixed=T)[[1]][1])
				
	}

	out = data.frame(time=time, year=year, month=month, day=day, hour=hour, min=min, sec=sec, lat=lat, long=long)

	return(out)
}
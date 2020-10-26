# Script to read in gpx files from Philippines 2008 fieldwork

setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
source("readGPX.R")
options(digits=22)

################################################################
## Concatentate GPX track files together into one spreadsheet
################################################################

files1 = list.files(path="gps/", pattern="^track[1][[:print:]]+gpx")
files2 = list.files(path="gps/", pattern="^track[2][[:print:]]+gpx")
gpsnum = c(rep(1, length(files1)), rep(2, length(files2)))
files = c(files1, files2)
len = length(files)
data = vector("list", len)

for(i in 1:len){
	infile = readGPX(paste("gps/",files[i],sep=""))
	data[i] = list(infile)
}

latlong = data.frame(data[1])
latlong$gpsnum = gpsnum[1]
for(i in 2:len){
	new = data.frame(data[i])
	new$gpsnum = gpsnum[i]
	latlong = merge(latlong, , all.x=T, all.y=T)
}

## Fill in data for #18 Pasil, Santander (missing)
pasillen = 79
pasillat = seq(9.419124, 9.415922, length.out = pasillen)
pasillong = seq(123.348109, 123.342953, length.out = pasillen)
pasilhour = c(rep(14, 13), rep(15,60), rep(16,6)) - 8 # to correct for GMT
pasilmin = formatC(c(47:59, 00:59, 0:5), width=2, flag="0")
pasiltime = paste("2008-08-31T", pasilhour, ":", pasilmin, ":30Z", sep="")
pasil = data.frame(time= pasiltime, year = 2008, month = 8, day = 31, hour = pasilhour, 
	min = as.numeric(pasilmin), sec=30, lat = pasillat, long = pasillong, gpsnum = 1)
latlong = merge(latlong, pasil, all.x=T, all.y=T)

## Fill in data for #22 Casay, Argao (missing)
casay = read.csv("gps/Casay survey.csv")
casay = casay[casay$use==1,]
casay = subset(casay, select=c(time, year, month, day, hour, min, sec, lat, long))
casay$gpsnum = 1
latlong = merge(latlong, casay, all.x=T, all.y=T)

# Sort the records by time
permut = order(latlong$year, latlong$month, latlong$day, latlong$hour, latlong$min, latlong$sec)
latlong = latlong[permut,]

write.csv(latlong, file=paste("gps/track.concat", Sys.Date(), ".csv", sep=""))


################################
## Match times to locations
################################

surv = read.csv("GPSSurveys.csv")

# Sygnathids

sygs = read.csv("GPSSygnathids.csv")
sygs$lat = NA
sygs$long = NA

len = dim(sygs)[1]

for(i in 1:len){
	#Get date and time information
	survey = sygs$Survey[i]
	survindex = which(surv$SurveyNum == survey)
	date = as.character(surv$Date[survindex])
	datesplit = strsplit(date,"/", fixed=T)[[1]]
	month = as.numeric(datesplit[1])
	day = as.numeric(datesplit[2])
	time = as.character(sygs$Time[i])
	timesplit = strsplit(time, ":", fixed=T)[[1]]
	hour = as.numeric(timesplit[1])
	min = as.numeric(timesplit[2])

	# Convert time to GMT
	hour = hour - 8
	if(hour <0){
		day = day-1
		hour = hour + 24
	}

	# Find the location records that match the date/time stamp
	latlongindex = which(latlong$month == month & latlong$day == day & latlong$hour == hour & latlong$min == min)

	# Calculate the mean lat/long for this minute
	sygs$lat[i] = mean(latlong$lat[latlongindex])
	sygs$long[i] = mean(latlong$long[latlongindex])
}

write.csv(sygs, file=paste("GPSSygnathids_wlatlong", Sys.Date(), ".csv", sep=""))



######################################################
## Trim track files to surveys only and output as GPX
######################################################
source("writeGPXheader.R")
source("writeGPXtrk.R")

outfilename = paste("gps/track.surveys.", Sys.Date(), ".gpx", sep="")

surv = read.csv("GPSSurveys.csv")


len = dim(surv)[1]

writeGPXheader(outfilename)

for(i in 1:len){
	if(surv$LinearFishSurvey[i] ==1 | surv$MappingFishSurvey[i]==1){
		# Find the survey segment
		survnum = surv$SurveyNum[i]
		date = strsplit(as.character(surv$Date[i]),"/", fixed=T)[[1]]
		month = as.numeric(date[1])
		day = as.numeric(date[2])
		starttime = strsplit(as.character(surv$StartTime[i]), ":", fixed=T)[[1]]
		endtime = strsplit(as.character(surv$EndTime[i]), ":", fixed=T)[[1]]
		starthour = as.numeric(starttime[1])
		startmin = as.numeric(starttime[2])
		endhour = as.numeric(endtime[1])
		endmin = as.numeric(endtime[2])

		# Convert time to GMT
		starthour = starthour - 8
		endhour = endhour - 8
		if(starthour <0 & endhour <0){
			day = day-1
			starthour = starthour + 24
			endhour = endhour + 24
		}
		if(starthour < 0 & endhour > 0){
			warning("This survey crosses between 2 days in GMT!")
			print(paste("Survey #", i, sep=""))
			## I'd then have to adjust the code to deal with this.... hopefully not an issue
		}

		startline = which(latlong$month == month & latlong$day == day & latlong$hour == starthour & latlong$min == startmin)
		startline = head(startline, n=1)
		endline = which(latlong$month == month & latlong$day == day & latlong$hour == endhour & latlong$min == endmin)
		endline = tail(endline, n=1)

		if(length(startline)>0 & length(endline)>0){

			lat = latlong$lat[startline:endline]
			long = latlong$long[startline:endline]
			time = latlong$time[startline:endline]
			hour = latlong$hour[startline:endline]
			min = latlong$min[startline:endline]

			# Average the positions within each minute
			j = 1
			k = 1
			avelen = length(unique(paste(hour,min)))
			latave = numeric(avelen)
			longave = numeric(avelen)
			timeave = character(avelen)
			while(j <= length(lat)){
				thishour = hour[j]
				thismin = min[j]
				m = which(hour==thishour & min == thismin)
				latave[k] = mean(lat[m])
				longave[k] = mean(long[m])
				timeave[k] = sub("[[:digit:]]{2}Z" , "30Z", as.character(time[m[1]]))
				k = k+1
				j = tail(m, n=1) + 1
			}
		
	
			# Write out the survey segment in GPX
			writeGPXtrk(outfilename, survnum, latave, longave, timeave)
		} else {
			warning(paste("No track for Survey #", i, sep=""))
		}	
	}
}

cat("</gpx>", file=outfilename, append=T)


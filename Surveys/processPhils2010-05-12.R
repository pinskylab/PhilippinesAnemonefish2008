


################################################################
## Concatentate GPX track files together into one spreadsheet
################################################################
# Script to read in gpx files from Philippines 2008 fieldwork
setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
source("readGPX.R")
options(digits=22)


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
	latlong = merge(latlong, new, all.x=T, all.y=T)
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

## Fill in data for #20 Poblacion, Argao (battery died after 7:03GMT, survey till 7:17GMT)
i = tail(which(latlong$month == 9 & latlong$day==6 & latlong$hour == 7 & latlong$min <= 17), 1)
poblen = 17-3
poblat = tail(seq(latlong$lat[i], 9.869093, length.out=poblen+1), poblen)
poblong = tail(seq(latlong$long[i], 123.596608, length.out=poblen+1), poblen)
pobhour = rep(7, poblen) 
pobmin = formatC(4:17, width=2, flag="0")
pobtime = paste("2008-08-31T", pobhour, ":", pobmin, ":30Z", sep="")
pob = data.frame(time= pobtime, year = 2008, month = 9, day = 6, hour = pobhour, 
	min = as.numeric(pobmin), sec=30, lat = poblat, long = poblong, gpsnum = 2)
latlong = merge(latlong, pob, all.x=T, all.y=T)

## Fill in data for #22 Casay, Argao (missing)
casay = read.csv("gps/trackCasay.csv")
casay = casay[casay$use==1,]
casay = subset(casay, select=c(time, year, month, day, hour, min, sec, lat, long))
casay$gpsnum = 1
latlong = merge(latlong, casay, all.x=T, all.y=T)

# Sort the records by time
permut = order(latlong$year, latlong$month, latlong$day, latlong$hour, latlong$min, latlong$sec)
latlong = latlong[permut,]

write.csv(latlong, file=paste("gps/track.concat", Sys.Date(), ".csv", sep=""))


################################
## Collections: Match times to locations
################################

surv = read.csv("GPSSurveys.csv")


##### Sygnathids

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


#####  Tissue Samples
#setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Analysis")
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")

surv = read.csv("GPSSurveys2009-01-08.csv")
data = read.csv("GPSSurveys2008-11-10.data.csv", as.is = c("ObsTime", "Spp", "Size1", "ID1"))
latlong = read.csv("gps/track.concat2008-11-21.csv")


# Find the largest pair of fish on each anemone
data$Rank1 = NA # column of the largest individual (ID1 through ID6)
data$Rank2 = NA
k = which(!is.na(data$Size1))
for(i in k){
	if(data$Size6[i] != ""){
		temp = as.numeric(unlist(strsplit(as.character(data$Size6[i]), split=",", fixed=T)))
	} else {
		temp = NA
	}
	j = sort(c(data$Size1[i],data$Size2[i],data$Size3[i],data$Size4[i],data$Size5[i],temp), index.return=T, decreasing=T)$ix 
	data$Rank1[i] = j[1]
	data$Rank2[i] = j[2]
}


collections = data.frame(SurveyNum = character(0), ObsTime = character(0), Spp = character(0), Size1 = numeric(0), ID1 = character(0), Notes=character(0), TopTwo = logical(0))
names = names(collections)

# add all fish that were sampled (ID is not blank)
k = which(data$ID1 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size1, ID1, Notes))
x$TopTwo = data$Rank1[k] == 1 | data$Rank2[k] == 1
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID2 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size2, ID2, Notes))
x$TopTwo = data$Rank1[k] == 2 | data$Rank2[k] == 2
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID3 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size3, ID3, Notes))
x$TopTwo = data$Rank1[k] == 3 | data$Rank2[k] == 3
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID4 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size4, ID4, Notes))
x$TopTwo = data$Rank1[k] == 4 | data$Rank2[k] == 4
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID5 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size5, ID5, Notes))
x$TopTwo = data$Rank1[k] == 5 | data$Rank2[k] == 5
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID6 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size6, ID6, Notes))
x$TopTwo = data$Rank1[k] == 6 | data$Rank2[k] == 6
names(x)= names
collections = rbind(collections, x)
dim(collections)

names(collections) = c("SurveyNum", "Time", "Spp", "Size", "ID", "Notes", "TopTwo")
collections$lat = NA
collections$long = NA
collections$Notes = as.character(collections$Notes)

len = dim(collections)[1]
for(i in 1:len){
	#Get date and time information
	survey = collections$Survey[i]
	survindex = which(surv$SurveyNum == survey)
	date = as.character(surv$Date[survindex])
	datesplit = strsplit(date,"/", fixed=T)[[1]]
	month = as.numeric(datesplit[1])
	day = as.numeric(datesplit[2])
	time = as.character(collections$Time[i])
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
	collections$lat[i] = mean(latlong$lat[latlongindex])
	collections$long[i] = mean(latlong$long[latlongindex])
}

# fix IDs 1-6 and 7-8 where we didn't label tubes before diving (insert average lat,long, time of the collections)
i = collections$ID == "1 to 6"
lat = mean(collections$lat[i])
long = mean(collections$long[i])
time = format(round(mean(strptime(collections$Time[i], "%H:%M")), units="mins"), "%H:%M") # mean time by closest minute
collections$Time[i] = time
collections$lat[i] = lat
collections$long[i] = long
collections$ID[i] = seq(1,6)
collections$Size[i] = NA

i = collections$ID == "7 to 8"
lat = mean(collections$lat[i])
long = mean(collections$long[i])
time = format(round(mean(strptime(collections$Time[i], "%H:%M")), units="mins"), "%H:%M") # mean time by closest minute
collections$Time[i] = time
collections$lat[i] = lat
collections$long[i] = long
collections$ID[i] = seq(7,8)
collections$Size[i] = NA

# fix Survey #14 where the GPS battery died partway through: insert a lat/long from a little after the survey
i = latlong$month==8 & latlong$day == 30 & latlong$hour == 5 & latlong$min == 31
lat = mean(latlong$lat[i], na.rm=T)
long = mean(latlong$long[i], na.rm=T)
j= collections$ID == 71 | collections$ID == 72
collections$lat[j] = lat
collections$long[j] = long
collections$Notes[j] = "GPS battery died. Latlong was recorded a short time after the sample was taken"


# fix Survey #22 (Casay): samples taken after the GPS track ended: (insert ave lat/long)
i = collections$SurveyNum == 22 & is.na(collections$lat)
collections$lat[i] = 9.832977396
collections$long[i] = 123.5688171
collections$Notes[i] = "Samples taken after survey ended. Use ave lat/long of this circular survey"


# fix survey #38: GPS was off until transect ended: use position from 9/17/08 at 06:52GMT
i = latlong$month==9 & latlong$day==17 & latlong$hour ==6 & latlong$min == 52
j = collections$SurveyNum == 38
collections$lat[j] = mean(latlong$lat[i], na.rm=T)
collections$long[j] = mean(latlong$long[i], na.rm=T)
collections$Notes[j] = "GPS was off. Use position recorded after survey ended"

# should also carry over notes field!

collections$Size = as.numeric(collections$Size)
permut = order(collections$SurveyNum, collections$Time)
collections = collections[permut,]

write.csv(collections, file=paste("Collections", Sys.Date(), ".csv", sep=""))



# Summarize collections by site (for BFAR)
# Add site, municipality and province
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Analysis")
collections = read.csv("Collections2008-11-21.csv")
surv = read.csv("GPSSurveys2008-11-10.csv")

collections = merge(collections, subset(surv, select=c("SurveyNum", "SiteNum", "Municipality", "Province")))
sitenums = unique(collections$SiteNum)
len = length(sitenums)
collsbysite = data.frame(SiteNum = sitenums, Municipality = character(len), Province = character(len), NumPRBI = numeric(len), NumAPCL = numeric(len), stringsAsFactors=FALSE)

for(i in 1:len){
	k = which(collections$SiteNum == sitenums[i])
	collsbysite$NumAPCL[i] = sum(collections$Spp[k] == "APCL")
	collsbysite$NumPRBI[i] = sum(collections$Spp[k] == "PRBI")
	collsbysite$Municipality[i] = as.character(unique(collections$Municipality[k]))
	collsbysite$Province[i]=as.character(unique(collections$Province[k]))
}

permut = order(collsbysite$SiteNum)
collsbysite = collsbysite[permut,]

write.csv(collsbysite, file = paste("Collections", Sys.Date(), ".bysite.csv", sep=""))




################################################################
## Trim track files to surveys only, clean, and output as GPX
################################################################
source("writeGPXheader.R")
source("writeGPXtrk.R")
source("writeGPXtrkseg.R")

outfilename = paste("gps/track.surveys.", Sys.Date(), ".gpx", sep="")

latlong = read.csv("gps/track.concat2008-10-24.csv")
surv = read.csv("GPSSurveys2008-10-24.csv")
cleaning = read.csv("Survey_digressions2008-10-24.csv")


len = dim(surv)[1]

writeGPXheader(outfilename)

for(i in 1:len){
	cat(paste(i, "\n", sep=""))
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
		thisgps = surv$GPSNum[i]
		pause = surv$Discontinuous[i]

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

		if(pause == 1){ # Discontinuous survey
			pausestarttime = strsplit(as.character(surv$PauseStart[i]), ":", fixed=T)[[1]]
			pauseendtime = strsplit(as.character(surv$PauseEnd[i]), ":", fixed=T)[[1]]
			pausestarthour = as.numeric(pausestarttime[1])
			pausestartmin = as.numeric(pausestarttime[2])
			pauseendhour = as.numeric(pauseendtime[1])
			pauseendmin = as.numeric(pauseendtime[2])

			# Convert time to GMT
			pausestarthour = pausestarthour - 8
			pauseendhour = pauseendhour - 8
			if(pausestarthour <0 & pauseendhour <0){
				pausestarthour = pausestarthour + 24
				pauseendhour = pauseendhour + 24
			}

			pausestartline = which(latlong$month == month & latlong$day == day & latlong$hour == pausestarthour & latlong$min == pausestartmin)
			pausestartline = tail(pausestartline , n=1)
			pauseendline = which(latlong$month == month & latlong$day == day & latlong$hour == pauseendhour & latlong$min == pauseendmin)
			pauseendline = head(pauseendline , n=1)		
		}

		if(length(startline)>0 & length(endline)>0){ # if we found start/end times for this survey in our latlong data

			if(pause==0){
				lat = latlong$lat[startline:endline]
				long = latlong$long[startline:endline]
				time = latlong$time[startline:endline]
				hour = latlong$hour[startline:endline]
				min = latlong$min[startline:endline]
				gpsnum = latlong$gpsnum[startline:endline]
			}	

			if(pause==1){ # if this is a discontinuous survey
				lat = latlong$lat[startline:pausestartline]
				long = latlong$long[startline:pausestartline]
				time = latlong$time[startline:pausestartline]
				hour = latlong$hour[startline:pausestartline]
				min = latlong$min[startline:pausestartline]
				gpsnum = latlong$gpsnum[startline:pausestartline]

				lat2 = latlong$lat[pauseendline:endline]
				long2 = latlong$long[pauseendline:endline]
				time2 = latlong$time[pauseendline:endline]
				hour2 = latlong$hour[pauseendline:endline]
				min2 = latlong$min[pauseendline:endline]
				gpsnum2 = latlong$gpsnum[pauseendline:endline]

				m = which(gpsnum2==thisgps)
				lat2 = lat2[m]
				long2 = long2[m]
				time2 = time2[m]
				hour2 = hour2[m]
				min2 = min2[m]
			}

			# Trim to the correct GPS
			m = which(gpsnum==thisgps)
			lat = lat[m]
			long = long[m]
			time = time[m]
			hour = hour[m]
			min = min[m]

			# Clean up digressions from linear in the linear surveys, based on the times in cleaning
			digr = survnum == cleaning$surveynum
			numdigr = sum(digr)
			if(numdigr>0){
				digrindex = which(digr)
				for(j in 1:numdigr){
					digrstart = strsplit(as.character(cleaning$start[digrindex[j]]), ":", fixed=T)[[1]]
					digrend = strsplit(as.character(cleaning$end[digrindex[j]]), ":", fixed=T)[[1]]
					digrstarthour = as.numeric(digrstart[1])
					digrstartmin = as.numeric(digrstart[2])
					digrendhour = as.numeric(digrend[1])
					digrendmin = as.numeric(digrend[2])

					# (Time is already in GMT in the cleaning data.frame)

					digrstartindex = which(hour == digrstarthour & min == digrstartmin)
					digrstartindex = tail(digrstartindex , n=1)
					digrendindex = which(hour == digrendhour & min == digrendmin)
					digrendindex = head(digrendindex , n=1)

					seg1 = TRUE # flag to process segment one
					if(pause == 1 & length(digrstartindex)==0){ # if not found in segment #1
						digrstartindex = which(hour2 == digrstarthour & min2 == digrstartmin)
						digrstartindex = tail(digrstartindex , n=1)
						digrendindex = which(hour2 == digrendhour & min2 == digrendmin)
						digrendindex = head(digrstartindex , n=1)

						if(length(digrstartindex)>0 & length(digrendindex)>0){	
							len = length(lat2)
							lat2 = lat2[c(1:digrstartindex, digrendindex:len)]
							long2 = long2[c(1:digrstartindex, digrendindex:len)]
							time2 = time2[c(1:digrstartindex, digrendindex:len)]
							hour2 = hour2[c(1:digrstartindex, digrendindex:len)]
							min2 = min2[c(1:digrstartindex, digrendindex:len)]
						}
						seg1 = FALSE # don't process segment one, since we just processed segment 2
					} 	
					if(length(digrstartindex)>0 & length(digrendindex)>0 & seg1 == T){	
						len = length(lat)
						lat = lat[c(1:digrstartindex, digrendindex:len)]
						long = long[c(1:digrstartindex, digrendindex:len)]
						time = time[c(1:digrstartindex, digrendindex:len)]
						hour = hour[c(1:digrstartindex, digrendindex:len)]
						min = min[c(1:digrstartindex, digrendindex:len)]
					}
					if(length(digrstartindex)==0 & length(digrendindex)==0){
						warning("No start/end time found for this digression")
					}
				}
			}

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

			if(pause == 1){
				j = 1
				k = 1
				avelen = length(unique(paste(hour2,min2)))
				latave2 = numeric(avelen)
				longave2 = numeric(avelen)
				timeave2 = character(avelen)
				while(j <= length(lat2)){
					thishour = hour2[j]
					thismin = min2[j]
					m = which(hour2==thishour & min2 == thismin)
					latave2[k] = mean(lat2[m])
					longave2[k] = mean(long2[m])
					timeave2[k] = sub("[[:digit:]]{2}Z" , "30Z", as.character(time2[m[1]]))
					k = k+1
					j = tail(m, n=1) + 1
				}
			}
	

	
			# Write out the survey segment in GPX
			if(pause==0){
				writeGPXtrk(outfilename, survnum, latave, longave, timeave)
			}
			if(pause==1){
				writeGPXtrk(outfilename, paste(survnum,".1",sep=""), latave, longave, timeave)
				writeGPXtrk(outfilename, paste(survnum,".2",sep=""), latave2, longave2, timeave2)
			}
			
		} else {
			warning(paste("No track for Survey #", i, sep=""))
		}
	
	}
}

cat("</gpx>", file=outfilename, append=T)



#########################################################################
###### Calculate density by survey using area calculations from ArcGIS
###### Calc Adult density for all spp
#########################################################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
#setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)


data = read.csv("GPSSurveys2008-11-10.data.csv")
surv = read.csv("GPSSurveys2009-01-08.csv")
area = read.csv("surveys.2010-05-12.area.csv")

surv = merge(surv, subset(area, select=c(SurveyNum, length, area, lat, long)), by="SurveyNum", all.x=T)

# anemone density, fish density, anemone area, fish biomass

# Add field names
anemlist = c("HECR", "ENQD", "STME", "MADO", "HEAR", "HEMG", "STHD", "STGI")
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")
for(thisspp in c(anemlist, fishlist)){
	field = paste("count", thisspp, sep="")
	surv[[field]] = NA
	field = paste("dens", thisspp, sep="")
	surv[[field]] = NA
	
	if(thisspp=="APCL"){
		field = paste("count", thisspp, "ad", sep="")
		surv[[field]] = NA
		field = paste("dens", thisspp, "ad", sep="")
		surv[[field]] = NA
		field = paste("count", thisspp, "juv", sep="")
		surv[[field]] = NA
		field = paste("dens", thisspp, "juv", sep="")
		surv[[field]] = NA
	}
}


# Find the largest pair on each anemone
data$Rank1 = 0 # column of the largest individual (ID1 through ID6)
data$Rank2 = 0
k = which(!is.na(data$Size1)) # makes our search shorter
for(i in k){
	if(data$Size6[i] != ""){
		temp = as.numeric(unlist(strsplit(as.character(data$Size6[i]), split=",", fixed=T)))
	} else {
		temp = NA
	}
	j = sort(c(data$Size1[i],data$Size2[i],data$Size3[i],data$Size4[i],data$Size5[i],temp), index.return=T, decreasing=T)$ix 
	data$Rank1[i] = j[1]
	data$Rank2[i] = j[2]
	if(is.na(data$Rank2[i])){ data$Rank2[i] = 0} # get rid of NAs
}

# Summarize count and density
for(thissurvey in surv$SurveyNum){
	for(thisspp in anemlist){
		count = sum(data$SurveyNum == thissurvey & data$AnemSpp == thisspp)
		k = surv$SurveyNum == thissurvey
		field = paste("count", thisspp, sep="")
		surv[[field]][k] = count		
		field = paste("dens", thisspp, sep="")
		surv[[field]][k] = count/surv$area[k]	
	}

	for(thisspp in fishlist){
		# from Spp fields
		k = (data$SurveyNum == thissurvey & data$Spp == thisspp)
		sizes = c(data$Size1[k], data$Size2[k], data$Size3[k], data$Size4[k], data$Size5[k]) 
		sizes2 = as.numeric(unlist(strsplit(as.character(data$Size6[k]), ",")))
		count = sum(!is.na(sizes)) + length(sizes2) # total number of fish		

		# from Spp2 fields
		k = (data$SurveyNum == thissurvey & data$Spp2 == thisspp)
		sizes3 = data$Spp2Size1[k]
		sizes4 = as.numeric(unlist(strsplit(as.character(data$Spp2Size2[k]), ",")))
		count3 = sum(!is.na(sizes3)) # from Spp2Size1 field
		count4 = length(sizes4) # from Spp2Size2 field

		k = surv$SurveyNum == thissurvey
		field = paste("count", thisspp, sep="")
		surv[[field]][k] = sum(count, count3, count4)		
		field = paste("dens", thisspp, sep="")
		surv[[field]][k] = sum(count, count3, count4)/surv$area[k]

		# adult and juv for APCL
		if(thisspp == "APCL"){
			# add adult dens and count for APCL
			k = (data$SurveyNum == thissurvey & data$Spp == thisspp)
			rank1 = c(data$Rank1[k]==1, data$Rank1[k]==2,data$Rank1[k]==3,data$Rank1[k]==4,data$Rank1[k]==5)
			rank2 = c(data$Rank2[k]==1, data$Rank2[k]==2,data$Rank2[k]==3,data$Rank2[k]==4,data$Rank2[k] == 5)
			ad = sum((rank1 | rank2) & (sizes >= 8)) # number of adults
			if(6 %in% rank1 | 6 %in% rank2){ # if one of the adult pair is in the Size6 column
				print(paste("Found adult in Size6 col: thispp=", thisspp, "thissurvey=", thissurvey, sep=""))
			}
			k = surv$SurveyNum == thissurvey
			field = paste("count", thisspp, "ad", sep="")
			surv[[field]][k] = ad		
			field = paste("dens", thisspp, "ad", sep="")
			surv[[field]][k] = ad/surv$area[k]
		
			# add juv dens and count for APCL
			countjuv = sum(c(sizes, sizes2, sizes3, sizes4)<=6, na.rm=T) # total number of fish
			k = surv$SurveyNum == thissurvey
			field = paste("count", thisspp, "juv", sep="")
			surv[[field]][k] = countjuv	
			field = paste("dens", thisspp, "juv", sep="")
			surv[[field]][k] = countjuv/surv$area[k]
		}
	}
}

write.csv(surv, paste("surveys", Sys.Date(),".density.csv", sep=""))





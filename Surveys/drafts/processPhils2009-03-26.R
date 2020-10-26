# Script to read in gpx files from Philippines 2008 fieldwork



################################################################
## Concatentate GPX track files together into one spreadsheet
################################################################
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
###### Calc Total and Adult density for APCL
#########################################################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
#setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)


data = read.csv("GPSSurveys2008-11-10.data.csv")
surv = read.csv("GPSSurveys2009-01-08.csv")
area = read.csv("surveys.2008-10-25.area.csv")

surv = merge(surv, subset(area, select=c(SurveyNum, length, area, lat)), by="SurveyNum", all.x=T)

# anemone density, fish density, anemone area, fish biomass

# Add field names
anemlist = c("HECR", "ENQD", "STME", "MADO", "HEAR", "HEMG", "STHD", "STGI")
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")
for(thisspp in c(anemlist, fishlist)){
	field = paste("count", thisspp, sep="")
	surv[[field]] = NA
	field = paste("dens", thisspp, sep="")
	surv[[field]] = NA
}
surv$countAPCLad = NA # adults and juveniles for APCL
surv$densAPCLad = NA

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
		
		# add adult dens and count for APCL
		if(thisspp=="APCL"){
			rank1 = c(data$Rank1[k]==1, data$Rank1[k]==2,data$Rank1[k]==3,data$Rank1[k]==4,data$Rank1[k]==5)
			rank2 = c(data$Rank2[k]==1, data$Rank2[k]==2,data$Rank2[k]==3,data$Rank2[k]==4,data$Rank2[k] == 5)
			ad = sum((rank1 | rank2) & (sizes >= 8)) # number of adults
			if(6 %in% rank1 | 6 %in% rank2){ # if one of the adult pair is in the Size6 column
				print(paste("Found adult in Size6 col: thispp=", thisspp, "thissurvey=", thissurvey, sep=""))
			}
			k = surv$SurveyNum == thissurvey
			surv$countAPCLad[k] = ad		
			surv$densAPCLad[k] = ad/surv$area[k]
		}

		# from Spp2 fields
		k = (data$SurveyNum == thissurvey & data$Spp2 == thisspp)
		count3 = sum(!is.na(data$Spp2Size1[k])) # from Spp2Size1 field
		count4 = length(unlist(strsplit(as.character(data$Spp2Size2[k]), ","))) # from Spp2Size2 field

		k = surv$SurveyNum == thissurvey
		field = paste("count", thisspp, sep="")
		surv[[field]][k] = sum(count, count3, count4)		
		field = paste("dens", thisspp, sep="")
		surv[[field]][k] = sum(count, count3, count4)/surv$area[k]
	}
}

write.csv(surv, paste("surveys", Sys.Date(),".density.csv", sep=""))


##############################################
## Extrapolate density to population size ####
##############################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")

surv = read.csv("surveys2009-03-31.density.csv")

reefarea = data.frame(name = c("Philippines", "Study Area", "CebuLeyte", "Cebu", "Leyte"), ArcArea = c(1260, 79,47, 18.9, 28.5), PolArea = c(5238, 498, 121, 85, 36), PolPerimArea = c(2418, 182, 35.4, 27.5, 7.95), TotalAreaPolArea = c(6499, 576, 168, 104, 64), TotalAreaPolPerim = c(3713, 264, 83, 47, 37))  # numbers are from PhilsReefArea 090318.R with reef arcs 150m wide

# survey numbers that are associated with either reef polygons or reef arcs in ReefBase GIS
polsurveys = c(2,3,4,5,6,7,8,40,41,42,43,48,73,74,75,76,78,94,95,96,97,98,99,100,101,102) 
arcsurveys = c(15,18,21,23,24,28,29,37,39,51,52,53,54,68,71,77,88,90,92,93)

# Total fish: using all surveys
cbind(reefarea$name, round(reefarea[,c("TotalAreaPolArea","TotalAreaPolPerim")]*mean(surv$densAPCL[surv$LinearFishSurvey==1 | surv$MappingFishSurvey==1], na.rm=T)*1000^2))   # convert fish/m2 to fish/km2

# Total fish: matching surveys to arcs and polygons
polsurvs = match(polsurveys, surv$SurveyNum)
arcsurvs = match(arcsurveys, surv$SurveyNum)
NumAPCL_polarea = round((reefarea$PolArea*mean(surv$densAPCL[polsurvs], na.rm=T) + reefarea$ArcArea*mean(surv$densAPCL[arcsurvs], na.rm=T))*1000^2) # using whole polygon areas
NumAPCL_polperim = round((reefarea$PolPerimArea*mean(surv$densAPCL[polsurvs], na.rm=T) + reefarea$ArcArea*mean(surv$densAPCL[arcsurvs], na.rm=T))*1000^2) # using only polygon perimeters (150m)

format(data.frame(reefarea$name,NumAPCL_polarea, NumAPCL_polperim), big.mark=",")
c("CebuLeyte per km", NumAPCL_polarea[3]/600, NumAPCL_polperim[3]/600)
c("Cebu per km", NumAPCL_polarea[4]/600, NumAPCL_polperim[4]/600)
c("Leyte per km", NumAPCL_polarea[5]/600, NumAPCL_polperim[5]/600)

# Adults: matching surveys to arcs and polygons and converting back to adults/km with 600 km coastline
polsurvs = match(polsurveys, surv$SurveyNum)
arcsurvs = match(arcsurveys, surv$SurveyNum)
NumAPCL_polarea = round((reefarea$PolArea*mean(surv$densAPCLad[polsurvs], na.rm=T) + reefarea$ArcArea*mean(surv$densAPCLad[arcsurvs], na.rm=T))*1000^2) # using whole polygon areas
NumAPCL_polperim = round((reefarea$PolPerimArea*mean(surv$densAPCLad[polsurvs], na.rm=T) + reefarea$ArcArea*mean(surv$densAPCLad[arcsurvs], na.rm=T))*1000^2) # using only polygon perimeters (150m)
out = data.frame(name=reefarea$name,NumAPCL_polarea, NumAPCL_polperim)
format(out, big.mark=",", digits=5)
c("CebuLeyte per km", NumAPCL_polarea[3]/600, NumAPCL_polperim[3]/600)
c("Cebu per km", NumAPCL_polarea[4]/600, NumAPCL_polperim[4]/600)
c("Leyte per km", NumAPCL_polarea[5]/600, NumAPCL_polperim[5]/600)

# bootstrap estimate of se of fish/km
#p = reefarea$PolArea[3]*surv$densAPCLad[polsurvs]*1000^2 # Area: total adults in polygons estimates (vary w/ density estimates)
#a = reefarea$ArcArea[3]*surv$densAPCLad[arcsurvs]*1000^2
#p = reefarea$PolPerimArea[3]*surv$densAPCLad[polsurvs]*1000^2 # Perim: total adults in polygons estimates 
#a = reefarea$ArcArea[3]*surv$densAPCLad[arcsurvs]*1000^2
#km=600
p = reefarea$PolPerimArea[4]*surv$densAPCLad[polsurvs]*1000^2 # Perim: Cebu
a = reefarea$ArcArea[4]*surv$densAPCLad[arcsurvs]*1000^2
km=250
#p = reefarea$PolPerimArea[5]*surv$densAPCLad[polsurvs]*1000^2 # Perim: Leyte 
#a = reefarea$ArcArea[5]*surv$densAPCLad[arcsurvs]*1000^2
#km=350
np = length(p)
na = length(a)
b = numeric(0)
c = numeric(0)
for(i in 1:10000){
	b = c(b,(mean(sample(p,np, replace=T))+mean(sample(a,na, replace=T)))/km) # adults/km
	c = c(c,(mean(sample(p,np, replace=T))+mean(sample(a,na, replace=T)))) # adults
}
mean(b)
sd(b, na.rm=T) # bootstrap estimate of se (resampling takes care of converting sd to se)
quantile(b, probs=c(0.025, 0.975))
mean(c)
sd(c, na.rm=T) # bootstrap estimate of se (resampling takes care of converting sd to se)
quantile(c, probs=c(0.025, 0.975))

# Adults/km directly from random site surveys
# ave
mean(surv$densAPCLad[surv$RandomSite==1], na.rm=T)*150*1000
sd(surv$densAPCLad[surv$RandomSite==1]*150*1000, na.rm=T)/sqrt(sum(!is.na(surv$densAPCLad[surv$RandomSite==1])))

# cebu
i = surv$RandomSite==1 & surv$Region =="Cebu"
mean(surv$densAPCLad[i], na.rm=T)*150*1000
sd(surv$densAPCLad[i]*150*1000, na.rm=T)/sqrt(sum(!is.na(surv$densAPCLad[i])))

# leyte
i = surv$RandomSite==1 & surv$Region =="Leyte"
mean(surv$densAPCLad[i], na.rm=T)*150*1000
sd(surv$densAPCLad[i]*150*1000, na.rm=T)/sqrt(sum(!is.na(surv$densAPCLad[i])))

i = surv$RandomSite==1 & surv$Region =="Cebu"
j = surv$RandomSite==1 & surv$Region =="Leyte"
t.test(surv$densAPCLad[i], surv$densAPCLac[j])



#########################################################################
###### Examine APCL sizes
#########################################################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
#setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")

data = read.csv("GPSSurveys2008-11-10.data.csv")
surv = read.csv("GPSSurveys2009-01-08.csv")
data = merge(data, subset(surv, select=c(SurveyNum, Region)))

k = (data$Spp == "APCL" & data$Region=="Cebu")
sizes_cebu = c(data$Size1[k], data$Size2[k], data$Size3[k], data$Size4[k], data$Size5[k])
sizes_cebu = c(sizes_cebu,as.numeric(unlist(strsplit(as.character(data$Size6[k]), ",")))) # from Size6 field
k = (data$Spp2 == "APCL" & data$Region=="Cebu")
sizes_cebu = c(sizes_cebu,data$Spp2Size1[k]) # from Spp2Size1 field
sizes_cebu = c(sizes_cebu,as.numeric(unlist(strsplit(as.character(data$Spp2Size2[k]), ",")))) # from Spp2Size2 field
sizes_cebu = sizes_cebu[!is.na(sizes_cebu)]

k = (data$Spp == "APCL" & data$Region=="Leyte")
sizes_leyte = c(data$Size1[k], data$Size2[k], data$Size3[k], data$Size4[k], data$Size5[k])
sizes_leyte = c(sizes_leyte,as.numeric(unlist(strsplit(as.character(data$Size6[k]), ",")))) # from Size6 field
k = (data$Spp2 == "APCL" & data$Region=="Leyte")
sizes_leyte = c(sizes_leyte,data$Spp2Size1[k]) # from Spp2Size1 field
sizes_leyte = c(sizes_leyte,as.numeric(unlist(strsplit(as.character(data$Spp2Size2[k]), ",")))) # from Spp2Size2 field
sizes_leyte = sizes_leyte[!is.na(sizes_leyte)]

quartz(height=8, width=5)
par(mfrow=c(3,1))
xlim = c(0, max(c(sizes_cebu, sizes_leyte)))
hist(c(sizes_cebu, sizes_leyte), main="All APCL", xlab="Length (cm)", xlim=xlim, breaks=15)
hist(sizes_cebu, main="Cebu APCL", xlab="Length (cm)", xlim=xlim, breaks=15)
abline(v = mean(sizes_cebu), lty=2)
hist(sizes_leyte, main="Leyte APCL", xlab="Length (cm)", xlim=xlim, breaks=15)
abline(v = mean(sizes_leyte), lty=2)

t.test(sizes_cebu, sizes_leyte)

# look for clusters
library(mclust)

clust = Mclust(sizes)
clust
clust$BIC
clust$parameters
par(mfcol=c(2,2))
plot(clust, sizes)

clust = mclustBIC(sizes)
sum = summary(clust, sizes)
sum
sum$parameters
par(mfrow=c(2,2))
plot(clust)
mclust1Dplot(sizes, classification = sum$classification, parameters=sum$parameters, what="classification")
mclust1Dplot(sizes, classification = sum$classification, parameters=sum$parameters, what="density")
abline(v = sum$parameters$mean, lty=3)
mclust1Dplot(sizes, classification = sum$classification, parameters=sum$parameters, what="uncertainty", uncertainty=sum$uncertainty)


# Check the collected samples for breeders (top 2 on anemone) that are too small (Hattori & Yanagisawa 1991 paper)
collected = read.csv("Collections2009-03-26.csv")
k = collected$Spp=="APCL"
summary(collected$Size[k & collected$TopTwo])
sum(collected$Size[k]<7 & collected$TopTwo[k]) # top 2 on anomene, but <7cm (not breeders)
sum(collected$Size[k]>7 & collected$TopTwo[k]) # top 2 on anomene, and >7cm (Real Breeders)
sum(collected$Size[k]<7 & !collected$TopTwo[k]) # not top 2 on anomene, but <7cm (real non-breeders)
sum(collected$Size[k]>7 & !collected$TopTwo[k]) # not top 2 on anomene, and >7cm (big non-breeders)


#############
### Plots
#############

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
#setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)

surv=read.csv("surveys2009-03-31.density.csv")

####### Density vs. lat, vis, depth, region, anemones #######
# Random Sites: PRBI and APCL
par(mfrow = c(2,4))
# PRBI
k = surv$RandomSite == 1
#par(mfrow = c(2,2))
plot(surv$lat[k], surv$densPRBI[k])
plot(surv$Visibility[k], surv$densPRBI[k])
plot((surv$DepthTop[k]+ surv$DepthBottom[k])/2, surv$densPRBI[k])
plot(surv$Region[k], surv$densPRBI[k])

require(vioplot)
k = !is.na(surv$densPRBI) & surv$RandomSite == 1
#vioplot(surv$densPRBI[surv$Region=="Cebu" & k], surv$densPRBI[surv$Region=="Danajon" & k], 
#	surv$densPRBI[surv$Region=="Leyte" & k], col="grey", names=c("Cebu", "Danajon", "Leyte"))

# APCL
k = surv$RandomSite == 1
#par(mfrow = c(2,2))
plot(surv$lat[k], surv$densAPCL[k])
plot(surv$Visibility[k], surv$densAPCL[k])
plot((surv$DepthTop[k]+ surv$DepthBottom[k])/2, surv$densAPCL[k])
plot(surv$Region[k], surv$densAPCL[k])

require(vioplot)
k = !is.na(surv$densAPCL) & surv$RandomSite == 1
vioplot(surv$densAPCL[surv$Region=="Cebu" & k], surv$densAPCL[surv$Region=="Danajon" & k], 
	surv$densAPCL[surv$Region=="Leyte" & k], col="grey", names=c("Cebu", "Danajon", "Leyte"))



# APCL: dist from city, anemones

# A. clarkii:  regress fish/m2 against citydist and anems (RANDOM sites)
quartz(height=4, width=8)
par(mfrow=c(1,3))
k = surv$Region == "Cebu" & surv$RandomSite==1
surv$citydist[k] = abs(surv$lat[k] - (10+17/60))
k = surv$Region == "Leyte" & surv$RandomSite==1
surv$citydist[k] =  abs(surv$lat[k] - 11.005)
k = (surv$Region == "Cebu" | surv$Region == "Leyte") & surv$RandomSite ==1
plot(surv$citydist[k], surv$densAPCL[k], xlab = "Distance from City (° Lat)", ylab = "Fish/m2", main="Distance from City")
m = lm(surv$densAPCL[k] ~ surv$citydist[k])
summary(m)
j = sort(surv$citydist[k], index.return=T)$ix
lines(surv$citydist[k][j], m$fitted.values[j])

k = (surv$Region == "Cebu" | surv$Region == "Leyte") & surv$RandomSite ==1
surv$densANEM[k] = surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k]
plot(surv$densANEM[k], surv$densAPCL[k], xlab = "Anemone/m2", ylab = "Fish/m2", main="Anemones")
m = lm(surv$densAPCL[k] ~ surv$densANEM[k])
summary(m)
j = sort(surv$densANEM[k], index.return=T)$ix
lines(surv$densANEM[k][j], m$fitted.values[j])

k = (surv$Region == "Cebu" | surv$Region == "Leyte") & surv$RandomSite ==1
plot(surv$citydist[k], surv$densANEM[k], xlab = "Distance from City (° Lat)", ylab = "Anemones/m2", main="Anemones and Distance")
m = lm(surv$densANEM[k] ~ surv$citydist[k])
summary(m)
j = sort(surv$citydist[k], index.return=T)$ix
lines(surv$citydist[k][j], m$fitted.values[j])


k = (surv$Region == "Cebu" | surv$Region == "Leyte") & surv$RandomSite ==1
m1 = lm(surv$densAPCL[k] ~ surv$citydist[k] + surv$densANEM[k])
m2 = lm(surv$densAPCL[k] ~ surv$citydist[k]) 
anova(m1, m2)
summary(m1)

m3 = lm(surv$densAPCL[k] ~ surv$citydist[k]*surv$densANEM[k])
summary(m3)
anova(m3, m2)


####### APCL histograms
par(mfrow=c(1,2))
k = surv$RandomSite == 1
breaks = seq(0, 0.030, by=0.005)
hist(surv$densAPCL[k & surv$Region == "Cebu"], breaks=breaks)
hist(surv$densAPCL[k & surv$Region == "Leyte"], breaks=breaks)



######## Plots of density vs. sitenum ##########

## Fish/m2 for APCL, PRBI
quartz(width=9.5,height=7)
par(mfrow=c(2,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu"
plot(surv$SiteNum[k], surv$densAPCL[k], main = "A. clarkii in Cebu", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.1))
k = surv$Region == "Leyte"
plot(surv$SiteNum[k], surv$densAPCL[k], main = "A. clarkii in Leyte", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.1))
k = surv$Region == "Danajon"
plot(surv$SiteNum[k], surv$densAPCL[k], main = "A. clarkii on Danajon Bank", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.1))

k = surv$Region == "Cebu"
plot(surv$SiteNum[k], surv$densPRBI[k], main = "P. biaculeatus in Cebu", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.03))
k = surv$Region == "Leyte"
plot(surv$SiteNum[k], surv$densPRBI[k], main = "P. biaculeatus in Leyte", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.03))
k = surv$Region == "Danajon"
plot(surv$SiteNum[k], surv$densPRBI[k], main = "P. biaculeatus on Danajon Bank", xlab = "Site Number", ylab = "fish per square meter", ylim=c(0,0.03))

## Fish/anem for APCL, PRBI
quartz(width=9.5,height=7)
par(mfrow=c(2,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu"
y = surv$densAPCL[k]/(surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k])
plot(surv$SiteNum[k], y, main = "A. clarkii in Cebu", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,4))
k = surv$Region == "Leyte"
y = surv$densAPCL[k]/(surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k])
plot(surv$SiteNum[k], y, main = "A. clarkii in Leyte", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,4))
k = surv$Region == "Danajon"
y = surv$densAPCL[k]/(surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k])
plot(surv$SiteNum[k], y, main = "A. clarkii in Danajon Bank", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,4))

k = surv$Region == "Cebu"
y = surv$densPRBI[k]/surv$densENQD[k]
plot(surv$SiteNum[k], y, main = "P. biaculeatus in Cebu", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,2))
k = surv$Region == "Leyte"
y = surv$densPRBI[k]/surv$densENQD[k]
plot(surv$SiteNum[k], y, main = "P. biaculeatus in Leyte", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,2))
k = surv$Region == "Danajon"
y = surv$densPRBI[k]/surv$densENQD[k]
plot(surv$SiteNum[k], y, main = "P. biaculeatus in Danajon Bank", xlab = "Site Number", ylab = "fish per anemone", ylim=c(0,2))


######## A. clarkii Plots of density vs. latitude ##########
printMeanSE <- function(x){ print(paste("Mean:", round(mean(x, na.rm=T)), "    SE:", round(sd(x, na.rm=T)/sqrt(sum(!is.na(x))), digits=3), sep="")) }

printMeanCI <- function(x){
	mu = mean(x, na.rm=T)
	se = sd(x, na.rm=T)/sqrt(sum(!is.na(x)))
	u = mu+1.96*se
	l = mu-1.96*se
	print(paste("Mean:", round(mu), "    CI:", round(l, digits=2), "-", round(u, digits=2), sep=""))
}

# a resampling approach. remember, though: means are normally distributed, even if data are not
printMeanCIboot <- function(x){
	mu = mean(x, na.rm=T)
	n = length(x)
	b = numeric(0)
	for(i in 1:10000){
		b = c(b, mean(sample(x, n, replace=T), na.rm=T))
	}	
	ci = quantile(b, probs = c(0.025, 0.975))
	print(paste("Mean:", round(mu), "    CI:", round(ci[1], digits=2), "-", round(ci[2], digits=2), sep=""))
}

# Fish/m2 for A clarkii, by latitude
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu"
plot(surv$lat[k], surv$densAPCL[k], main = "A. clarkii in Cebu", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.1))
k = surv$Region == "Leyte"
plot(surv$lat[k], surv$densAPCL[k], main = "A. clarkii in Leyte", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.1))
k = surv$Region == "Danajon"
plot(surv$lat[k], surv$densAPCL[k], main = "A. clarkii on Danajon Bank", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.1))

# Fish/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 14000)
k = surv$Region == "Cebu"
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Cebu\n(150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Leyte"
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Leyte\n(150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Danajon"
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii on Danajon Bank\n(150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)

# Random sites: Fish/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 4700)
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Cebu\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Leyte\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii on Danajon Bank\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)

# All BUT random sites: Fish/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 14000)
k = surv$Region == "Cebu" & surv$RandomSite != 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Cebu\n(all but random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Leyte" & surv$RandomSite != 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii in Leyte\n(all but random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)
k = surv$Region == "Danajon" & surv$RandomSite != 1
plot(surv$lat[k], surv$densAPCL[k]*150*1000, main = "A. clarkii on Danajon Bank\n(all but random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densAPCL[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCL[k]*150*1000)


# Adults/m2 for A clarkii, by latitude
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu"
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults in Cebu", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.03))
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)
k = surv$Region == "Leyte"
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults in Leyte", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.03))
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)
k = surv$Region == "Danajon"
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults on Danajon Bank", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=c(0,0.03))
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)

# Random sites: Adults/m2 for A clarkii, by latitude
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims=c(0,max(surv$densAPCLad[surv$RandomSite==1]))
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults in Cebu", xlab = "Latitude (°)", ylab = "fish per square meter", ylim=ylims)
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)
k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults in Leyte", xlab = "Latitude (°)", ylab = "fish per square meter", ylim= ylims)
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)
k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k], main = "A. clarkii adults on Danajon Bank", xlab = "Latitude (°)", ylab = "fish per square meter", ylim= ylims)
mean(surv$densAPCLad[k], na.rm=T); sd(surv$densAPCLad[k], na.rm=T)


# Adults/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu"
ylims = c(0, 4000)
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Cebu\n(150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)
k = surv$Region == "Leyte"
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Leyte\n(150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)
k = surv$Region == "Danajon"
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults on Danajon Bank\n(150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)

# Random sites: Adults/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 1200)
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Cebu\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)

k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Leyte\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)

k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults on Danajon Bank\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)


# Random sites: Adults/km for A clarkii, by latitude (150m wide reef): Cebu & Leyte together
k = (surv$Region == "Leyte" | surv$Region == "Cebu") & surv$RandomSite == 1
printMeanSE(surv$densAPCLad[k]*150*1000)
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults on Cebu and Leyte\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)


# All BUT random sites: Adults/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
k = surv$Region == "Cebu" & surv$RandomSite != 1
ylims = c(0, 4000)
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Cebu\n(all but random, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)
k = surv$Region == "Leyte" & surv$RandomSite != 1
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults in Leyte\n(all but random, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)
k = surv$Region == "Danajon" & surv$RandomSite != 1
plot(surv$lat[k], surv$densAPCLad[k]*150*1000, main = "A. clarkii adults on Danajon Bank\n(all but random, 150m reef)", xlab = "Latitude (°)", ylab = "adults per km", ylim= ylims)
abline(h = mean(surv$densAPCLad[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densAPCLad[k]*150*1000)


# FOR EVOL PAPER: A clarkii fish and adults/m2 for Cebu and Leyte, Random and Non-Random
a = surv[order(surv$lat),]
quartz(width=7, height=3.5)
par(mfrow=c(1,2))
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
ylims = c(0, 0.1)
xlims = c(9.4, 11.4)
k = a$Region == "Cebu" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Cebu", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])
sd(a$densAPCL[k])/sum(k)
k = a$Region == "Cebu" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)

xlims = c(9.95, 11.1)
k = a$Region == "Leyte" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Leyte", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])
sd(a$densAPCL[k])/sum(k)
k = a$Region == "Leyte" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)





################
### Plot A. clarkii against distance from city

# A. clarkii: regress density against distance (deg lat) from cebu city (all sites)
quartz(height=4, width=7)
par(mfrow=c(1,2))
k = surv$Region == "Cebu" & !is.na(surv$lat)
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main ="Cebu")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & !is.na(surv$lat)
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main="Leyte")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

legend("topright", legend=c("South", "North"), pch=c(16,1))


# A. clarkii:  regress density against distance (deg lat) from cebu city (RANDOM sites)
quartz(height=4, width=7)
par(mfrow=c(1,2))
k = surv$Region == "Cebu" & !is.na(surv$lat) & surv$RandomSite==1
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main ="Cebu")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & !is.na(surv$lat) & surv$RandomSite==1
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main="Leyte")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

legend("topright", legend=c("South", "North"), pch=c(16,1))



# A. clarkii:  regress density against distance (deg lat) from cebu or ormoc city (all sites)
quartz(height=4, width=7)
par(mfrow=c(1,2))
k = surv$Region == "Cebu" & !is.na(surv$lat)
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main ="Cebu")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & !is.na(surv$lat)
x = abs(surv$lat[k] - 11.005)
plot(x, surv$densAPCL[k], xlab = "Distance from Ormoc City (° Latitude)", ylab = "Fish/m2", main="Leyte")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

legend("topleft", legend=c("South", "North"), pch=c(16,1))


# A. clarkii:  regress density against distance (deg lat) from cebu or ormoc city (RANDOM sites)
quartz(height=4, width=7)
par(mfrow=c(1,2))
k = surv$Region == "Cebu" & !is.na(surv$lat) & surv$RandomSite==1
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densAPCL[k], xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/m2", main ="Cebu")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & !is.na(surv$lat) & surv$RandomSite==1
x = abs(surv$lat[k] - 11.005)
plot(x, surv$densAPCL[k], xlab = "Distance from Ormoc City (° Latitude)", ylab = "Fish/m2", main="Leyte")
i = surv$lat[k] < (11.005) # southern points
points(x[i], surv$densAPCL[k][i], pch=16)
m = lm(surv$densAPCL[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

legend("topleft", legend=c("South", "North"), pch=c(16,1))


# A. clarkii:  regress fish/anem against distance (deg lat) from cebu or ormoc city (RANDOM sites)
quartz(height=4, width=7)
par(mfrow=c(1,2))
k = surv$Region == "Cebu" & !is.na(surv$lat) & surv$RandomSite==1
y = surv$countAPCL[k]/(surv$countHECR[k]+surv$countENQD[k]+surv$countHEAR[k]+surv$countHEMG[k]+surv$countMADO[k]+surv$countSTGI[k]+surv$countSTHD[k]+surv$countSTME[k])
x = abs(surv$lat[k] - (10+17/60))
plot(x, y, xlab = "Distance from Cebu City (° Latitude)", ylab = "Fish/anemone", main ="Cebu")
i = surv$lat[k] < (10+17/60) # southern points
points(x[i], y[i], pch=16)
m = lm(y ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & !is.na(surv$lat) & surv$RandomSite==1
y = surv$countAPCL[k]/(surv$countHECR[k]+surv$countENQD[k]+surv$countHEAR[k]+surv$countHEMG[k]+surv$countMADO[k]+surv$countSTGI[k]+surv$countSTHD[k]+surv$countSTME[k])
y[is.nan(y)]=0
x = abs(surv$lat[k] - 11.005)
plot(x, y, xlab = "Distance from Ormoc City (° Latitude)", ylab = "Fish/anemone", main="Leyte")
i = surv$lat[k] < (11.005) # southern points
points(x[i], y[i], pch=16)
m = lm(y ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

legend("topleft", legend=c("South", "North"), pch=c(16,1))

# A. clarkii: one figure: regress fish/anem against distance (deg lat) from cebu or ormoc city (RANDOM sites)
quartz(height=4, width=4)
par(mfrow=c(1,1))
k = surv$Region == "Cebu" & !is.na(surv$lat) & surv$RandomSite==1
y = surv$countAPCL[k]/(surv$countHECR[k]+surv$countENQD[k]+surv$countHEAR[k]+surv$countHEMG[k]+surv$countMADO[k]+surv$countSTGI[k]+surv$countSTHD[k]+surv$countSTME[k])
x = abs(surv$lat[k] - (10+17/60))
k = surv$Region == "Leyte" & !is.na(surv$lat) & surv$RandomSite==1
y = c(y, surv$countAPCL[k]/(surv$countHECR[k]+surv$countENQD[k]+surv$countHEAR[k]+surv$countHEMG[k]+surv$countMADO[k]+surv$countSTGI[k]+surv$countSTHD[k]+surv$countSTME[k]))
y[is.nan(y)]=0
x = c(x, abs(surv$lat[k] - 11.005))
plot(x, y, xlab = "Distance from City (° Latitude)", ylab = "Fish/anemone", main="Fish/anemone vs. Distance from City")
m = lm(y ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])





######## P. biaculeatus Plots of density vs. latitude ##########

# Random sites: Fish/km for PRBI, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 400)
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Cebu\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)
k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Leyte\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)
k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus on Danajon Bank\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)

# All BUT Random sites: Fish/km for PRBI, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 5000)
k = surv$Region == "Cebu" & surv$RandomSite != 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Cebu\n(non-random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)
k = surv$Region == "Leyte" & surv$RandomSite != 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus in Leyte\n(non-random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)
k = surv$Region == "Danajon" & surv$RandomSite != 1
plot(surv$lat[k], surv$densPRBI[k]*150*1000, main = "P. biaculeatus on Danajon Bank\n(non-random sites, 150m reef)", xlab = "Latitude (°)", ylab = "fish per km", ylim= ylims)
abline(h = mean(surv$densPRBI[k]*150*1000, na.rm=T), lty=3)
printMeanSE(surv$densPRBI[k]*150*1000)


######### A. clarkii boxplots by region #############

# Random sites: Adults/km for A clarkii, assuming 150m wide reef
#quartz(width=3, height=3)
quartz(width=8, height=8)
par(cex=2, cex.axis = 1.6, cex.lab = 1.7, bty="l", omi=c(0,0.1,0,0), lwd=2)

#boxplot(I(150*1000*densAPCLad) ~ Region, data=surv, subset = surv$RandomSite ==1, main = "A. clarkii adults\n(random sites, 150m reef)", ylab = "adults per km", range=0)

x = surv[surv$Region != "Danajon", ]
x$Region = x$Region[,drop=TRUE]
#title = "A. clarkii adults\n(random sites, 150m reef)" # for notebook
title = "Visual census density" # for ppt
boxplot(I(150*1000*densAPCLad) ~ Region, data=x, subset = x$RandomSite ==1, main = title, ylab = "adults per km", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=3, col="dark grey")

# Ave across Cebu & Leyte
quartz(width=7, height=8)
par(cex=2, cex.axis = 1.6, cex.lab = 1.7, bty="l", omi=c(0,0.1,0,0), lwd=2)

x = surv[surv$Region != "Danajon", ]
x$Region = " "
x$Region = x$Region[,drop=TRUE]
#title = "A. clarkii adults\n(random sites, 150m reef)" # for notebook
title = "Visual census density" # for ppt
boxplot(I(150*1000*densAPCLad)~Region, data=x, subset = x$RandomSite ==1, main = title, ylab = "adults per km", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=3, col="dark grey")


######################################################################
###### Summarize fish-anemone associations
######################################################################

# fish/anemone by spp (ave and stdev)
# number of anems seen with each spp of fish (and %)
# number of fish seen with each spp of anem (and %)

setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)
source("superpose.eb.R")

data = read.csv("GPSSurveys2008-11-07.data.csv") # survey data

anemlist = c("HECR", "ENQD", "STME", "MADO", "HEAR", "HEMG", "STHD", "STGI")
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")


# Summarize fish/anem by fish spp (conditional on fish present on anemone)
symbiosis = matrix(nrow = length(anemlist), ncol = 4*length(fishlist))
rownames(symbiosis) = anemlist
colnames(symbiosis) = paste(rep(fishlist, rep(4, length(fishlist))), c(".numanem", ".numfish", ".mean", ".sd"), sep="")

for(thisspp in anemlist){
	for(thisfish in fishlist){
		k = which(data$AnemSpp == thisspp & data$Spp == thisfish)
		if(length(k)>0){
			count = numeric(0)
			for(i in k){
				count1 = sum(!is.na(c(data$Size1[i], data$Size2[i], data$Size3[i], data$Size4[i], data$Size5[i])))
				count2 = length(unlist(strsplit(as.character(data$Size6[i]), ","))) # from Size6 field
				count3 = sum(!is.na(data$Spp2Size1[i])) # from Spp2Size1 field
				count4 = length(unlist(strsplit(as.character(data$Spp2Size2[i]), ","))) # from Spp2Size2 field
				count=c(count,count1+count2+count3+count4)
			}
			symbiosis[which(thisspp==anemlist), 4*which(thisfish==fishlist)-3] = length(count)
			symbiosis[which(thisspp==anemlist), 4*which(thisfish==fishlist)-2] = sum(count)
			symbiosis[which(thisspp==anemlist), 4*which(thisfish==fishlist)-1] = mean(count)
			symbiosis[which(thisspp==anemlist), 4*which(thisfish==fishlist)] = sd(count)
		}
	}
}

# make a plot of ave fish/anem, by fish spp.
windows(11,6)
par(mfrow=c(2,4))
for(i in 1:length(fishlist)){
	x.abcis = barplot(symbiosis[,4*i-1], col="black", ylim=c(0,8), ylab="Fish/anem", xlab="Anemone species", 
		main=fishlist[i], cex.names=0.5)
	superpose.eb(x.abcis, symbiosis[,4*i-1], symbiosis[,4*i], lwd=2, col="grey")
}


# plot of anem use by fish and plot of fish hosting by anems
windows(11,8)
par(mfrow=c(2,1))

y = symbiosis[,seq(2,dim(symbiosis)[2], by=4)]
denom = colSums(y, na.rm=T)
for(i in 1:length(denom)){
	y[,i] = y[,i]/denom[i]*100
}
colnames(y) = fishlist
barplot(y, beside=T, col=seq(0,length(anemlist)-1), ylim=c(0,100), ylab="Probability of occurring on anemone (%)", xlab="Fish species", 
		main="Use of anemones by fish", cex.names=1)
legend("topright", legend=anemlist, fill=seq(0,length(anemlist)-1), ncol=3, cex=0.8)

y = t(symbiosis[,seq(1,dim(symbiosis)[2], by=4)])
denom = colSums(y, na.rm=T)
for(i in 1:length(denom)){
	y[,i] = y[,i]/denom[i]*100
}
rownames(y) = fishlist
barplot(y, beside=T, col=seq(0,length(fishlist)-1), ylim=c(0,100), ylab="Probability of hosting fish (%)", xlab="Anemone species", 
		main="Hosting of fish by anemones", cex.names=1)
legend("topright", legend=fishlist, fill=seq(0,length(fishlist)-1), ncol=1, cex=0.8)






########################
## Statistics

k = surv$RandomSite == 1
a = surv$densAPCL[k & surv$Region == "Cebu"]
b = surv$densAPCL[k & surv$Region == "Leyte"]
summary(a)
summary(b)
t.test(a,b)

# exclude Tangkaan Beach
k = surv$RandomSite == 1 & surv$SurveyNum != 59
a = surv$densAPCL[k & surv$Region == "Cebu"]
b = surv$densAPCL[k & surv$Region == "Leyte"]
summary(a)
summary(b)
t.test(a,b)

# all sites
a = surv$densAPCL[surv$Region == "Cebu"]
b = surv$densAPCL[surv$Region == "Leyte"]
summary(a)
summary(b)
t.test(a,b)

# all sites except Tangkaan
k = surv$SurveyNum != 59
a = surv$densAPCL[k & surv$Region == "Cebu"]
b = surv$densAPCL[k & surv$Region == "Leyte"]
summary(a)
summary(b)
t.test(a,b)

## Adults (Random Sites)
a = surv$densAPCLad[surv$RandomSite==1 & surv$Region == "Cebu"]*150*1000
b = surv$densAPCLad[surv$RandomSite==1 & surv$Region == "Leyte"]*150*1000
summary(a)
summary(b)
t.test(a,b)

## Adults (all but Random Sites)
a = surv$densAPCLad[surv$RandomSite!=1 & surv$Region == "Cebu"]*150*1000
b = surv$densAPCLad[surv$RandomSite!=1 & surv$Region == "Leyte"]*150*1000
summary(a)
summary(b)
t.test(a,b)
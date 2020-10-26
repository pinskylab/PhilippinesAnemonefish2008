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


## Tissue Samples
setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")

surv = read.csv("GPSSurveys2008-11-07.csv")
data = read.csv("GPSSurveys2008-11-07.data.csv", as.is = c("ObsTime", "Spp", "Size1", "ID1"))
latlong = read.csv("gps/track.concat2008-10-22.csv")

collections = data.frame(SurveyNum = character(0), ObsTime = character(0), Spp = character(0), Size1 = numeric(0), ID1 = character(0))
names = names(collections)

k = which(data$ID1 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size1, ID1))
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID2 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size2, ID2))
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID3 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size3, ID3))
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID4 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size4, ID4))
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID5 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size5, ID5))
names(x)= names
collections = rbind(collections, x)
dim(collections)

k = which(data$ID6 !="")
length(k)
x<-subset(data[k,], select=c(SurveyNum, ObsTime, Spp, Size6, ID6))
names(x)= names
collections = rbind(collections, x)
dim(collections)

collections$lat = NA
collections$long = NA

len = dim(collections)[1]
for(i in 1:len){
	#Get date and time information
	survey = collections$Survey[i]
	survindex = which(surv$SurveyNum == survey)
	date = as.character(surv$Date[survindex])
	datesplit = strsplit(date,"/", fixed=T)[[1]]
	month = as.numeric(datesplit[1])
	day = as.numeric(datesplit[2])
	time = as.character(collections$ObsTime[i])
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



write.csv(collections, file=paste("Collections", Sys.Date(), ".csv", sep=""))

# Summarize collections by site (for BFAR)
# Add site, municipality and province

write.csv(colcsbysite, file = paste("Collections", Sys.Date(), ".bysite.csv", sep=""))




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
#########################################################################
setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)


data = read.csv("GPSSurveys2008-10-26.data.csv")
surv = read.csv("GPSSurveys2008-10-26.csv")
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
		k = (data$SurveyNum == thissurvey & data$Spp == thisspp)
		count = sum(!is.na(c(data$Size1[k], data$Size2[k], data$Size3[k], data$Size4[k], data$Size5[k])))
		count2 = length(unlist(strsplit(as.character(data$Size6[k]), ","))) # from Size6 field
		count3 = sum(!is.na(data$Spp2Size1[k])) # from Spp2Size1 field
		count4 = length(unlist(strsplit(as.character(data$Spp2Size2[k]), ","))) # from Spp2Size2 field

		k = surv$SurveyNum == thissurvey
		field = paste("count", thisspp, sep="")
		surv[[field]][k] = sum(count, count2, count3, count4)		
		field = paste("dens", thisspp, sep="")
		surv[[field]][k] = sum(count, count2, count3, count4)/surv$area[k]	
	}
}

write.csv(surv, "surveys2008-10-26.density.csv")



#############
### Plots
#############

setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)

surv=read.csv("surveys2008-10-26.density.csv")


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


# APCL histograms
par(mfrow=c(1,2))
k = surv$RandomSite == 1
breaks = seq(0, 0.030, by=0.005)
hist(surv$densAPCL[k & surv$Region == "Cebu"], breaks=breaks)
hist(surv$densAPCL[k & surv$Region == "Leyte"], breaks=breaks)


# Crude map of random sites (APCL/PRBI)
par(mfrow=c(1,2))
expfact = 500
expstep = 0.5
k = surv$RandomSite == 1 & surv$Region == "Cebu"
cex = expfact*surv$densAPCL[k]+expstep
plot(rep(1,sum(k)),surv$lat[k], cex = cex, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APCL Density", ylab = "Latitude")
k = surv$RandomSite == 1 & surv$Region == "Leyte"
cex=expfact*surv$densAPCL[k]+expstep
pch = (surv$densAPCL[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$RandomSite == 1 & surv$Region == "Cebu"
cex = expfact*surv$densPRBI[k]+expstep
pch = (surv$densPRBI[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "PRBI Density", ylab = "Latitude")
k = surv$RandomSite == 1 & surv$Region == "Leyte"
cex=expfact*surv$densPRBI[k]+expstep
pch = (surv$densPRBI[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.001, 0.005, 0.01)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="fish/m2")

# Crude map of all sites (APCL/PRBI) in fish/m2
par(mfrow=c(1,2))
expfact = 500
expstep = 0.5
k = surv$Region == "Cebu"
cex = expfact*surv$densAPCL[k]+expstep
plot(rep(1,sum(k)),surv$lat[k], cex = cex, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APCL Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densAPCL[k]+expstep
pch = (surv$densAPCL[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
cex = expfact*surv$densPRBI[k]+expstep
pch = (surv$densPRBI[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "PRBI Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densPRBI[k]+expstep
pch = (surv$densPRBI[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.001, 0.005, 0.01)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="fish/m2")



# Crude map of all sites (APCL/PRBI) in fish/anemone
par(mfrow=c(1,2))
expfact = 5
expstep = 0.5
k = surv$Region == "Cebu"
y = surv$densAPCL[k]/(surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k])
cex = expfact*y+expstep
pch = (surv$densAPCL[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APCL Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densAPCL[k]/(surv$densHECR[k]+surv$densENQD[k]+surv$densHEAR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k]+surv$densSTHD[k]+surv$densSTME[k])
cex = expfact*y+expstep
pch = (surv$densAPCL[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
y = surv$densPRBI[k]/surv$densENQD[k]
cex = expfact*y+expstep
pch = (surv$densPRBI[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "PRBI Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densPRBI[k]/surv$densENQD[k]
cex = expfact*y+expstep
pch = (surv$densPRBI[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.1, 0.5, 1)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="fish/anemone")




# Crude map of all sites: APOC, APML, APSA, APPE
par(mfrow=c(2,2))
expfact = 500
expstep = 0.5
k = surv$Region == "Cebu"
cex = expfact*surv$densAPOC[k]+expstep
plot(rep(1,sum(k)),surv$lat[k], cex = cex, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APOC Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densAPOC[k]+expstep
pch = (surv$densAPOC[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
cex = expfact*surv$densAPML[k]+expstep
pch = (surv$densAPML[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APML Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densAPML[k]+expstep
pch = (surv$densAPML[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
cex = expfact*surv$densAPSA[k]+expstep
pch = (surv$densAPSA[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APSA Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densAPSA[k]+expstep
pch = (surv$densAPSA[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.001, 0.005, 0.01)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="fish/m2")

k = surv$Region == "Cebu"
cex = expfact*surv$densAPPE[k]+expstep
pch = (surv$densAPPE[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APPE Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densAPPE[k]+expstep
pch = (surv$densAPPE[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)


# Crude map of all sites: APOC, APML, APSA, APPE in fish/anemone
par(mfrow=c(2,2))
expfact = 5
expstep = 0.5
k = surv$Region == "Cebu"
y = surv$densAPOC[k]/(surv$densHEMG[k]+surv$densSTGI[k]+surv$densSTME[k])
cex = expfact*y+expstep
plot(rep(1,sum(k)),surv$lat[k], cex = cex, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APOC Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densAPOC[k]/(surv$densHEMG[k]+surv$densSTGI[k]+surv$densSTME[k])
cex = expfact*y+expstep
pch = (surv$densAPOC[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
y = surv$densAPML[k]/(surv$densENQD[k])
cex = expfact*y+expstep
pch = (surv$densAPML[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APML+APFR Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densAPML[k]/(surv$densENQD[k])
cex = expfact*y+expstep
pch = (surv$densAPML[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
y = surv$densAPSA[k]/(surv$densHECR[k]+surv$densSTME[k])
cex = expfact*y+expstep
pch = (surv$densAPSA[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APSA Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densAPSA[k]/(surv$densHECR[k]+surv$densSTME[k])
cex = expfact*y+expstep
pch = (surv$densAPSA[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.1, 0.5, 1)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="fish/anemone")

k = surv$Region == "Cebu"
y = surv$densAPPE[k]/(surv$densHECR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k])
cex = expfact*y+expstep
pch = (surv$densAPPE[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "APPE Density", ylab = "Latitude")
k = surv$Region == "Leyte"
y = surv$densAPPE[k]/(surv$densHECR[k]+surv$densHEMG[k]+surv$densMADO[k]+surv$densSTGI[k])
cex = expfact*y+expstep
pch = (surv$densAPPE[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)



# Crude map of all sites: HECR, ENQD, STME, HEMG
par(mfrow=c(2,2))
expfact = 500
expstep = 0.5
k = surv$Region == "Cebu"
cex = expfact*surv$densHECR[k]+expstep
plot(rep(1,sum(k)),surv$lat[k], cex = cex, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "HECR Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densHECR[k]+expstep
pch = (surv$densHECR[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
cex = expfact*surv$densENQD[k]+expstep
pch = (surv$densENQD[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "ENQD Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densENQD[k]+expstep
pch = (surv$densENQD[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

k = surv$Region == "Cebu"
cex = expfact*surv$densSTME[k]+expstep
pch = (surv$densSTME[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "STME Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densSTME[k]+expstep
pch = (surv$densSTME[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)

leg = c(0, 0.001, 0.005, 0.01)
legend("topright", legend=leg, pch=c(20,1,1,1,1), 
	pt.cex = expfact*leg+expstep, title="anemones/m2")

k = surv$Region == "Cebu"
cex = expfact*surv$densHEMG[k]+expstep
pch = (surv$densHEMG[k] == 0)*19+1
plot(rep(1,sum(k)),surv$lat[k], cex = cex, pch=pch, xlim=c(0.8, 2.2), ylim=c(9.2, 11.8), main = "HEMG Density", ylab = "Latitude")
k = surv$Region == "Leyte"
cex=expfact*surv$densHEMG[k]+expstep
pch = (surv$densHEMG[k] == 0)*19+1
points(rep(2,sum(k)),surv$lat[k], cex = cex, pch=pch)




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
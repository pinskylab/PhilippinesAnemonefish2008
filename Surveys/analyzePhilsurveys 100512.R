
####################################
## Summary Stats on the Surveys ####
####################################
setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
surv = read.csv("surveys2009-09-23.density.csv")

# distance
k = surv$RandomSite == 1
sum(surv$length[k])/1000 # total length in km
sum(surv$length[k])/1000/(252+223)*100 # % of Cebu+Leyte study area (from 1/8/10 calculation in ArcGIS)
mean(surv$length[k]) # mean distance of random surveys
sd(surv$length[k])/sqrt(sum(k))

# area
k = surv$RandomSite == 1
sum(surv$area[k])# total area in m2


# distance and area of reef surveys
	# survey numbers that are associated with either reef polygons or reef arcs in ReefBase GIS
polsurveys = c(2,3,4,5,6,7,8,40,41,42,43,48,73,74,75,76,78,94,95,96,97,98,99,100,101,102) 
arcsurveys = c(15,18,21,23,24,28,29,37,39,51,52,53,54,68,71,77,88,90,92,93)

polsurvs = match(polsurveys, surv$SurveyNum)
arcsurvs = match(arcsurveys, surv$SurveyNum)
k = c(polsurvs, arcsurvs)
sum(surv$length[k])/1000 # total length in km
sum(surv$length[k])/1000/(252+223)*100 # % of Cebu+Leyte study area (from 1/8/10 calculation in ArcGIS)

sum(surv$area[k])# total area in m2




# duration
dur = strptime(surv$EndTime, format="%H:%M") - strptime(surv$StartTime, format="%H:%M")
dur[surv$Discontinuous==1] = dur[surv$Discontinuous==1] - as.numeric(strptime(surv$PauseEnd[surv$Discontinuous==1], format="%H:%M") - strptime(surv$PauseStart[surv$Discontinuous==1], format="%H:%M")) # remove pauses
surv$duration = dur

k = surv$RandomSite == 1
mean(surv$duration[k]) # mean duration of random surveys
sd(surv$duration[k])/sqrt(sum(k))


# depth
k = surv$RandomSite == 1
par(mfrow=c(1,2))
hist(surv$DepthTop[k])
hist(surv$DepthBottom[k])

mean(surv$DepthTop[k])
sd(surv$DepthTop[k])/sqrt(sum(k))
mean(surv$DepthBottom[k])
sd(surv$DepthBottom[k])/sqrt(sum(k))

mean(surv$DepthBottom[k] - surv$DepthTop[k])



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


# What's the average size of our defined juvs (<= 6 cm) and adults (>= 8 cm)?
collected = read.csv("Collections2009-03-26.csv")
min(collected$Size, na.rm=T) # 2 cm min
max(collected$Size, na.rm=T) # 16 cm max
mean(collected$Size[collected$Size<=6], na.rm=T) # 5.06 cm juvs
sd(collected$Size[collected$Size<=6], na.rm=T)/sqrt(sum(collected$Size<=6, na.rm=T)) # +/- 0.078 cm juvs
mean(collected$Size[collected$Size>=8 & collected$TopTwo == TRUE], na.rm=T) # 10.03 cm adults
sd(collected$Size[collected$Size>=8 & collected$TopTwo == TRUE], na.rm=T)/sqrt(sum(collected$Size>=8 & collected$TopTwo == TRUE, na.rm=T)) # +/- 0.097 cm adults




#################################
### Plots
#################################

setwd("/Users/mpinsky/Documents/Stanford/Philippines/2008/Surveys")
#setwd("C:/Documents and Settings/Lab/Desktop/Users/Malin/Philippines/2008 Surveys")
options(digits=22)

surv=read.csv("surveys2010-05-12.density.csv")



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

# Fish/m2 for A clarkii, by latitude (all sites)
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

## Adult APCL

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


# Adults/km for A clarkii, by latitude, assuming 150m wide reef (all sites)
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

# NON-RANDOM FISH SURVEYS (random and mapping): Adults/m2 for A. clarkii by latitude
a = surv[order(surv$lat),]
quartz(width=7, height=3.5)
par(mfrow=c(1,2))
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
ylims = c(0, 0.1)

xlims = c(9.4, 11.4)
k = a$Region == "Cebu" & a$RandomSite != 1 & (a$LinearFishSurvey == 1 | a$MappingFishSurvey==1)
plot(a$lat[k], a$densAPCLad[k], main = "Cebu", xlab = "Latitude (°)", ylab = "adults/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # SE across non-random sites
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site

xlims = c(9.95, 11.1)
k = a$Region == "Leyte" & a$RandomSite != 1 & (a$LinearFishSurvey == 1 | a$MappingFishSurvey==1)
plot(a$lat[k], a$densAPCLad[k], main = "Leyte", xlab = "Latitude (°)", ylab = "adults/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # SE across non-random sites
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site


# NON-RANDOM LINEAR FISH SURVEYS: Adults/m2 for A. clarkii by latitude
a = surv[order(surv$lat),]
quartz(width=7, height=3.5)
par(mfrow=c(1,2))
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
ylims = c(0, 0.1)

xlims = c(9.4, 11.4)
k = a$Region == "Cebu" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCLad[k], main = "Cebu", xlab = "Latitude (°)", ylab = "adults/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # SE across non-random sites
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site

xlims = c(9.95, 11.1)
k = a$Region == "Leyte" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCLad[k], main = "Leyte", xlab = "Latitude (°)", ylab = "adults/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # SE across non-random sites
aggregate(a$densAPCLad[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site



## Juvenile APCL

# Random sites: Juvs/km for A clarkii, by latitude, assuming 150m wide reef
quartz(width=9.5, height=3)
par(mfrow=c(1,3))
par(cex=.7, cex.lab = 1.4)
ylims = c(0, 3000)
k = surv$Region == "Cebu" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLjuv[k]*150*1000, main = "A. clarkii juveniles in Cebu\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "juveniles per km", ylim= ylims)
abline(h = mean(surv$densAPCLjuv[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLjuv[k]*150*1000)

k = surv$Region == "Leyte" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLjuv[k]*150*1000, main = "A. clarkii juveniles in Leyte\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "juveniles per km", ylim= ylims)
abline(h = mean(surv$densAPCLjuv[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLjuv[k]*150*1000)

k = surv$Region == "Danajon" & surv$RandomSite == 1
plot(surv$lat[k], surv$densAPCLjuv[k]*150*1000, main = "A. clarkii juveniles in Danajon\n(random sites, 150m reef)", xlab = "Latitude (°)", ylab = "juveniles per km", ylim= ylims)
abline(h = mean(surv$densAPCLjuv[k]*150*1000), lty=3)
printMeanSE(surv$densAPCLjuv[k]*150*1000)






# FOR EVOL PAPER (SUBMISSION #1): A clarkii fish and adults/m2 for Cebu and Leyte, Random and Non-Random Linear Fish Surveys
a = surv[order(surv$lat),]
quartz(width=7, height=3.5)
par(mfrow=c(1,2))
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
ylims = c(0, 0.1)

xlims = c(9.4, 11.4)
k = a$Region == "Cebu" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Cebu", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # se across non-random sites
aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site
k = a$Region == "Cebu" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)

xlims = c(9.95, 11.1)
k = a$Region == "Leyte" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Leyte", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])  # mean density for non-random sites
sd(a$densAPCL[k])/sum(k) # SE
sum(k) # n
k = a$Region == "Leyte" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)


# FOR PHILS REPORT: A clarkii fish and adults/m2 for Danajon, Random and Non-Random
a = surv[order(surv$lat),]
quartz(width=10.5, height=6)
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
par(mfrow=c(1,3))
ylims = c(0, 0.1)

xlims = c(9.4, 11.4)
k = a$Region == "Cebu" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Cebu (east coast)", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])
sd(a$densAPCL[k])/sum(k)
k = a$Region == "Cebu" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)

xlims = c(9.95, 11.1)
k = a$Region == "Leyte" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Leyte (west coast)", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])
sd(a$densAPCL[k])/sum(k)
k = a$Region == "Leyte" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)

xlims = c(10.17, 10.3)
k = a$Region == "Danajon" & a$RandomSite != 1 & (a$LinearFishSurvey == 1)
plot(a$lat[k], a$densAPCL[k], main = "Bohol (Danajon Bank)", xlab = "Latitude (°)", ylab = "fish/m2", ylim= ylims, xlim=xlims, type="p", pch=4)
mean(a$densAPCL[k])
sd(a$densAPCL[k])/sum(k)
k = a$Region == "Danajon" & a$RandomSite == 1
lines(a$lat[k], a$densAPCL[k], type="o", pch=20, lty=2)
lines(a$lat[k], a$densAPCLad[k], type="o", pch=20)



################
### Plot A. clarkii or anem against distance from city

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
legend("topleft", legend=c("South", "North"), pch=c(16,1))

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


# Anem:  regress anem/m2 against citydist (RANDOM sites)
quartz(height=5, width=8)
par(mfrow=c(1,2))

k = surv$Region == "Cebu" & surv$RandomSite ==1
x = abs(surv$lat[k] - (10+17/60))
plot(x, surv$densANEM[k], xlab = "Distance from Cebu City (° Lat)", ylab = "Anemones/m2", main="Cebu")
m = lm(surv$densANEM[k] ~ x)
summary(m)
j = sort(x, index.return=T)$ix
lines(x[j], m$fitted.values[j])

k = surv$Region == "Leyte" & surv$RandomSite ==1
x = abs(surv$lat[k] - 11.005)
plot(x, surv$densANEM[k], xlab = "Distance from Ormoc City (° Lat)", ylab = "Anemones/m2", main="Leyte")
m = lm(surv$densANEM[k] ~ x)
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

# Like EVOL PAPER APCL graph: PRBI fish and adults/m2 for Cebu, Bohol, Leyte, Random and Non-Random Linear Fish Surveys vs. Longitude
a = surv[order(surv$long),]
quartz(width=7, height=3.5)
par(cex=.7, cex.lab = 1.7, cex.axis=1.4, mai=c(0.7,0.7,0.5,0.2))
ylims = c(0, 3)
sites = c(7,8,9,10,11,1,2,19)

xlims = c(123.6, 125.3)
k = a$RandomSite != 1 & (a$LinearFishSurvey == 1) & a$SiteNum %in% sites # non-random sites
plot(a$long[k], a$densPRBI[k]*100, main = "", xlab = "Longitude (°E)", ylab = "Fish/100 m2", ylim= ylims, xlim=xlims, type="p", pch=4)
aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean) # mean density for each non-random site
mean(aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean)$x) # mean density across non-random sites
sd(aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), mean)$x)/sqrt(length(unique(a$SiteNum[k]))) # se across non-random sites
aggregate(a$densAPCL[k], by=list(a$SiteNum[k]), sd)$x/sqrt(table(a$SiteNum[k])) # SE for each non-random site
sum(k) # n
table(a$SiteNum[k]) # n per site
k = a$RandomSite == 1 & a$SiteNum %in% sites 
lines(a$long[k], a$densPRBI[k]*100, type="o", pch=20, lty=2)


########### All spp vs. latitude ##########


# Random sites: Adults/m2, by latitude (all spp)
quartz(width=9.5, height=7)
par(mfrow=c(2,4))
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")
for(thisspp in fishlist){
	dens = paste("dens", thisspp, "ad", sep="")
	x = surv[surv$Region != "Danajon", c("Region", "RandomSite", "lat", dens)]
	x$Region = x$Region[,drop=TRUE]
	names(x) = c("Region", "RandomSite", "lat", "dens")
	plot(dens ~ lat, data=x, subset = x$RandomSite == 1 & x$Region=="Cebu", main = thisspp, ylab = "adults/m2", col="black", pch=16)
	points(dens ~ lat, data=x, subset = x$RandomSite ==1 & x$Region=="Leyte", main = thisspp, ylab = "adults/m2", col="red", pch=16)
}
plot(0,0, xaxt="n", yaxt="n", col="white", bty="n", xlab="", ylab="")
legend("topright", legend=c("Cebu", "Leyte"), pch=16, col=c("black", "red"), cex=2, bty="n")


# All sites: Adults/m2, by latitude (all spp)
quartz(width=9.5, height=7)
par(mfrow=c(2,4))
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")
for(thisspp in fishlist){
	dens = paste("dens", thisspp, "ad", sep="")
	x = surv[surv$Region != "Danajon", c("Region", "RandomSite", "lat", dens)]
	x$Region = x$Region[,drop=TRUE]
	names(x) = c("Region", "RandomSite", "lat", "dens")
	plot(dens ~ lat, data=x, subset = x$Region=="Cebu", main = thisspp, ylab = "adults/m2", col="black", pch=c(4,16)[(dens>0)+1])
	points(dens ~ lat, data=x, subset = x$Region=="Leyte", main = thisspp, ylab = "adults/m2", col="red", pch=c(4,16)[(dens>0)+1])
}
plot(0,0, xaxt="n", yaxt="n", col="white", bty="n", xlab="", ylab="")
legend("topright", legend=c("Cebu", "Leyte"), pch=16, col=c("black", "red"), cex=2, bty="n")

# All sites: Count, by latitude (all spp)
quartz(width=9.5, height=7)
par(mfrow=c(2,4))
fishlist = c("PRBI", "APCL", "APOC", "APML", "APSA", "APPE", "APPY")
for(thisspp in fishlist){
	dens = paste("count", thisspp, sep="")
	x = surv[surv$Region != "Danajon", c("Region", "RandomSite", "lat", dens)]
	x$Region = x$Region[,drop=TRUE]
	names(x) = c("Region", "RandomSite", "lat", "dens")
	plot(dens ~ lat, data=x, subset = x$Region=="Cebu", main = thisspp, ylab = "# fish per dive", col="black", pch=c(4,16)[(dens>0)+1])
	points(dens ~ lat, data=x, subset = x$Region=="Leyte", main = thisspp, col="red", pch=c(4,16)[(dens>0)+1])
}
plot(0,0, xaxt="n", yaxt="n", col="white", bty="n", xlab="", ylab="")
legend("topright", legend=c("Cebu", "Leyte"), pch=16, col=c("black", "red"), cex=2, bty="n")



######### Boxplots by region #############

# APCL: Random sites: Adults/km for A clarkii, assuming 150m wide reef
#quartz(width=3, height=3)
quartz(width=8, height=8)
par(cex=2, cex.axis = 1.6, cex.lab = 1.7, bty="l", omi=c(0,0.1,0,0), lwd=2)

#boxplot(I(150*1000*densAPCLad) ~ Region, data=surv, subset = surv$RandomSite ==1, main = "A. clarkii adults\n(random sites, 150m reef)", ylab = "adults per km", range=0)

x = surv[surv$Region != "Danajon", ]
x$Region = x$Region[,drop=TRUE]
#title = "A. clarkii adults\n(random sites, 150m reef)" # for notebook
title = "Visual census density" # for ppt
boxplot(I(150*1000*densAPCLad) ~ Region, data=x, subset = x$RandomSite ==1, main = title, ylab = "adults per km", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=3, col="dark grey")

# APCL: Ave across Cebu & Leyte
quartz(width=7, height=8)
par(cex=2, cex.axis = 1.6, cex.lab = 1.7, bty="l", omi=c(0,0.1,0,0), lwd=2)

x = surv[surv$Region != "Danajon", ]
x$Region = " "
x$Region = x$Region[,drop=TRUE]
#title = "A. clarkii adults\n(random sites, 150m reef)" # for notebook
title = "Visual census density" # for ppt
boxplot(I(150*1000*densAPCLad)~Region, data=x, subset = x$RandomSite ==1, main = title, ylab = "adults per km", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=3, col="dark grey")

### APCL (Cebu-Leyte) vs. PRBI: only Random Sites, Adults (= all PRBI)
prbisites = c(7,8,9,10,11,1,2,19) # PRBI sites to use
x = surv[,c('SurveyNum', 'LinearFishSurvey', 'RandomSite', 'Region', 'SiteNum', 'densAPCLad', 'densPRBI')]
x = reshape(x, direction="long", varying=c('densAPCLad', 'densPRBI'), v.names='dens', times=c('APCL', 'PRBI'), idvar = 'SurveyNum', timevar='Species')
x = x[(x$Species == 'APCL' & x$Region != 'Danajon') | (x$Species == 'PRBI' & x$SiteNum %in% prbisites),] # trim out PRBI sites that are outside the xsect, and trim out APCL Danajon
x = x[x$RandomSite == 1, ]
x$Regspp[x$Species=='APCL'] = paste(x$Region[x$Species=='APCL'], x$Species[x$Species=='APCL'], sep="")
x$Regspp[x$Species=='PRBI'] = 'PRBI'
boxplot(I(100*dens)~Regspp, data=x, main = "", ylab = "Adults/100 m2", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=1, col="dark grey")

summary(lm(dens ~ Regspp, data=x)) # p = 0.1705

summary(lm(dens ~ Regspp, data=x, subset=x$Regspp == 'CebuAPCL' | x$Regspp == 'PRBI')) # p = 0.087 # CebuAPCL vs. PRBI
t.test(dens ~ Regspp, data=x, subset=x$Regspp == 'CebuAPCL' | x$Regspp == 'PRBI') # p = 0.069

### APCL vs. PRBI: only Random Sites, Adults (= all PRBI)
x = surv[,c('SurveyNum', 'LinearFishSurvey', 'RandomSite', 'Region', 'SiteNum', 'densAPCLad', 'densPRBI')]
x = reshape(x, direction="long", varying=c('densAPCLad', 'densPRBI'), v.names='dens', times=c('APCL', 'PRBI'), idvar = 'SurveyNum', timevar='Species')
x = x[(x$Species == 'APCL' & x$Region != 'Danajon') | (x$Species == 'PRBI' & x$SiteNum %in% prbisites),] # trim out PRBI sites that are outside the xsect, and trim out APCL Danajon
x = x[x$RandomSite == 1, ]
x$Regspp[x$Species=='APCL'] = paste(x$Region[x$Species=='APCL'], x$Species[x$Species=='APCL'], sep="")
x$Regspp[x$Species=='PRBI'] = 'PRBI'
quartz(width=5, height=5)
#pdf(paste("figures/APCL+PRBI adults per m2 boxplot ", Sys.Date(), ".pdf", sep=""), width=5, height=5)
boxplot(I(1000*dens)~Species, data=x, main = "APCL Cebu+Leyte, PRBI East-West", ylab = "Adults/1000 m2", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=1, col="dark grey", pars=list(bty='n', font.axis=3), names=c("A. clarkii", "P. biaculeatus"))
#dev.off()

summary(mod<-lm(dens ~ Species, data=x)) # p = 0.162, PRBI = 0.000589, APCL = 0.00168 (2.9x higher than PRBI)
t.test(dens ~ Species, data=x) # p = 0.0706, PRBI = 0.000589, APCL = 0.00168 (APCL 2.9x higher than PRBI)
	# PRBI = 88.4 +/- 37.8 SE fish/km (adults/km since all fish are breeding adults?)
wilcox.test(dens ~ Species, data=x) # p = 0.373, PRBI = 0.000589, APCL = 0.00168 (2.9x higher than PRBI)
library(exactRankTests); wilcox.exact(dens ~ Species, data=x) # p = 0.38

anova(mod)
unique(round(as.numeric(fitted(mod)), digits=6))

summary(mod<-lm(dens ~ Species, data=x, subset=!(x$Species=='PRBI' & x$SiteNum==19))) # p = 0.162, PRBI = 0.000589, APCL = 0.00168 (2.9x higher than PRBI) # remove Site 19 for PRBI
t.test(dens ~ Species, data=x, subset=!(x$Species=='PRBI' & x$SiteNum==19)) # p = 0.02, PRBI = 0.000351, APCL = 0.00168 (APCL 4.8x higher than PRBI)	# PRBI = 52.7 +/- 14.4 SE fish/km (adults/km)
unique(round(as.numeric(fitted(mod)), digits=6))

summary(mod<-lm(dens ~ Species, data=x, subset=x$Region=='Cebu')) # p = 0.119 # only Cebu


### APCL vs. PRBI: all Sites, Adults (= all PRBI)
x = surv[,c('SurveyNum', 'LinearFishSurvey', 'RandomSite', 'Region', 'SiteNum', 'densAPCLad', 'densPRBI')]
x = reshape(x, direction="long", varying=c('densAPCLad', 'densPRBI'), v.names='dens', times=c('APCL', 'PRBI'), idvar = 'SurveyNum', timevar='Species')
x = x[(x$Species == 'APCL' & x$Region != 'Danajon') | (x$Species == 'PRBI' & x$SiteNum %in% prbisites),] # trim out PRBI sites that are outside the xsect, and trim out APCL Danajon
x$Regspp[x$Species=='APCL'] = paste(x$Region[x$Species=='APCL'], x$Species[x$Species=='APCL'], sep="")
x$Regspp[x$Species=='PRBI'] = 'PRBI'
boxplot(I(1000*dens)~Species, data=x, main = "", ylab = "Adults/1000 m2", range=0, boxwex=0.7, staplewex=0.2, outwex = 0.5, lwd=1, col="dark grey", parts=list(bty='n'))

summary(mod<-lm(dens ~ Species, data=x)) # p = 0.162


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

y = symbiosis[,seq(2,dim(symbiosis)[2], by=4)] # only pull numfish columns
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


# plot of anem use by APCL and PRBI
y = symbiosis[,c('PRBI.numfish', 'APCL.numfish')] # only pull numfish columns
denom = colSums(y, na.rm=T)
for(i in 1:length(denom)){
	y[,i] = y[,i]/denom[i]*100
}
y[is.na(y)] =0
cols = c('black', 'red', 'yellow', 'white', 'orange', 'green', 'purple', 'grey')
quartz(width=8, height=4)
barplot(y, beside=T, ylim=c(0,100), col=cols, ylab="Fish (%)", xlab="Fish species", main="Use of anemones by fish", cex.names=1, names=c('P. biaculeatus', 'A. clarkii'))
legend("topright", legend=c('H. crispa', 'E. quadricolor', 'S. mertensi', 'M. doreensis', 'H. aurora', 'H. magnifica', 'S. haddoni', 'S. gigantea'), fill=cols, ncol=2, cex=0.8, bty='n')

# plot of anem filling by APCL and PRBI
y = t(symbiosis[,seq(1,dim(symbiosis)[2], by=4)])
for(i in 1:ncol(y)){
	y[,i] = y[,i]/sum(y[,i], na.rm=T)*100
}
y[is.na(y)] =0
y = rbind(y, colSums(y[3:7,]))
y = y[c(1,2,8),]
cols = c('black', 'red', 'yellow')
quartz(width=8, height=4)
barplot(y, beside=T, ylim=c(0,100), col=cols, ylab="Anemones (%)", xlab="Anemone species", main="Filling of anemones by fish", cex.names=1, names=anemlist)
legend("topleft", legend=c('P. biaculeatus', 'A. clarkii', 'Other Amphiprion'), fill=cols, ncol=1, cex=0.8, bty='n')




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
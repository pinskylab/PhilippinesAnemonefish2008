# Simulation of a two-part dispersal kernel: 60% retained, rest spread uniformly to get
# an sd of 30km
max=100
x = runif(400, 0, max)
x = c(x, rep(0, 600))
sd(x)

hist(x)

# Gaussian: Difference between axial standard dev, mean dispersal distance,and median
x = rnorm(1000, sd=5)
median(abs(x))
mean(abs(x))

# Laplace: Difference between axial standard dev, mean dispersal distance,and median
x=rexp(1000, rate=.29) * (2*round(runif(1000, 0,1))-1)
sd(x)
mean(abs(x))


#################################################################
#### Can we have 60% self-recruitment to 500m reef and sd of 7km?
	# Set up a graph
	quartz(width=6, height=4)
	par(mfrow=c(1,2), cex=1)
	col1 = "grey90" # the background color
	col2 = "maroon" # the local color
	col3 = "steelblue" # the outside color


## Continuous Reef, constant density across space (dispersal kernel = arrival kernel)
	# Gaussian
	len = 100000
	x = rnorm(len, sd=11)
	sum(abs(x)<0.25)/len # 8km: only 2.5% retained on 500m reef, only 2.5% of settlers are local origin
	sum(abs(x)<0.5)/len # 8km: only 5% retained on 1km reef
	# 7km: 2.8% retained on 500m reef, 5.7% retained on 1km reef
	# 10km: 2% retained on 500m reef, 4% on 1km reef
	# 11km: 1.8% retained on 500m, 3.6% on 1km
	# 12km: 1.6% retained on 500m reef, 3.3% retained on 1km reef
	quantile(abs(x), c(0.5, 0.75, 0.95))
	mean(abs(x)) # mean dispersal is 9.6 km for 12km stdev
		# 8.8 km for 11 km stdev
		# 8 km for 10km stdev

	# try pnorm instead of random numbers
	pnorm(.25, sd=11)-pnorm(-.25, sd=11) # 1.8% retained on 500 m reef with 11km sigma


	
	# Laplace (from Wikipedia generating random variables)
	len = 100000
	sd = 11
	u = runif(len, -1/2, 1/2)
	x = -sd/sqrt(2)*sign(u)*log(1-2*abs(u))
	#hist(x)
	sd(x)
	sum(abs(x)<0.25)/len # only 5% retained, only 5% of settlers are local origin
	# 11km: 3% retained on 500m reef

	# Bateman leptokurtic (see Wright 1969)

	# Graph x
	minreef = floor(min(x)/10)*10-0.5
	maxreef = ceiling(max(x)/10)*10+0.5
	breaks=seq(minreef, maxreef, 1)
	cols = rep(col2, length(breaks-1))
	cols[breaks< (-0.5)] = col3
	cols[breaks>=0.5] = col3
	hist(x, breaks , col=cols, border=cols, main="Continuous", xlab="Distance (km)", freq=F)


## Fragmented Reef, 500m-1km focal reef, all others > 10km
	# Gaussian
	len = 10000
	x = rnorm(len, sd=7)
	sum(abs(x)<0.25)/len # 2% of dispersers are retained on 500m reef
	sum(abs(x)<0.5)/len # 6% of dispersers are retained on 1km reef
	sum(abs(x)<0.25)/(sum(abs(x)>10)+sum(abs(x)<0.25)) # 12% of settlers are local origin on 500m reef
	sum(abs(x)<0.5)/(sum(abs(x)>10)+sum(abs(x)<0.5)) # 27% of settlers are local origin on 1km reef

	# Laplace (from Wikipedia generating random variables)
	len = 100000
	sd = 7
	u = runif(len, -1/2, 1/2)
	x = -sd/sqrt(2)*sign(u)*log(1-2*abs(u))
	#hist(x)
	sd(x)
	sum(abs(x)<0.25)/len # only 5% retained on 500m reef
	sum(abs(x)<0.5)/len # only 10% retained on 1km reef
	sum(abs(x)<0.25)/(sum(abs(x)>10)+sum(abs(x)<0.25)) # 30% of settlers are local origin
	sum(abs(x)<0.5)/(sum(abs(x)>10)+sum(abs(x)<0.5)) # 42% of settlers are local origin



## Fragmented Reef, 1km focal reef, others 1km every 10km
	# Gaussian
	len = 100000
	x = rnorm(len, sd=12)
	sum(abs(x)<0.5)/len # 6% of dispersers are retained on 1km reef
	mx = max(abs(x)) # what's the max distance traveled?
	numreefs = round((mx-0.5)/10) # number of patches we need
	outside = 0 # number of immigrants from outside
	dists = x[abs(x)<(1/2)] # distances traveled by successful dispersers
	for(i in 1:numreefs){ # add up the immigrants from each patch reef
		j = abs(x)>(i*10-0.5) & abs(x)<(i*10+0.5)
		outside = outside + sum(j)
		dists = c(dists, x[j])
	}
	sum(abs(x)<1/2)/(outside+sum(abs(x)<1/2)) # 50% of settlers are local origin
	# 12km: 33% are local origin
	sd(x)
	sd(dists) # stdev of dispersal distances is 8km (or 12km)
	quantile(abs(dists), c(0.5, 0.75, 0.95))

	# Laplace (from Wikipedia generating random variables)
	len = 100000
	sd = 7
	u = runif(len, -1/2, 1/2)
	x = -sd/sqrt(2)*sign(u)*log(1-2*abs(u))
	#hist(x)
	sd(x)
	sum(abs(x)<0.5)/len # only 10% retained on 1km reef
	mx = max(abs(x)) # what's the max distance traveled?
	numreefs = round((mx-0.5)/10) # number of patches we need
	outside = 0 # number of immigrants from outside
	for(i in 1:numreefs){ # add up the immigrants from each patch reef
		thispatch = sum(abs(x)>(i*10-0.5) & abs(x)<(i*10+0.5))
		print(thispatch)
		outside = outside + thispatch
	}
	sum(abs(x)<1/2)/(outside+sum(abs(x)<1/2)) # 75% of settlers are local origin

	# Graph x
	minreef = floor(min(x)/10)*10-0.5
	maxreef = ceiling(max(x)/10)*10+0.5
	breaks= seq(minreef, maxreef, 1)
	cols=rep(rep(c(col3, col1), c(1,9)), (length(breaks)-2)/10)
	cols[breaks== -0.5] = col2
	hist(x, breaks, col=cols, border=cols, main="Discontinuous", xlab="Distance (km)", freq=F)

## Fragmented Reef, 0.5km focal reef, others 0.5km every X km
	# Gaussian
	len = 100000
	reef = 1.0/2
	reefdist = 15
	x = rnorm(len, sd=11)
	sum(abs(x)<reef)/len # 1.6% of dispersers are retained on 1km reef
	mx = max(abs(x)) # what's the max distance traveled?
	numreefs = round((mx-reef)/reefdist) # number of patches we need
	outside = 0 # number of immigrants from outside
	dists = x[abs(x)<(reef)] # distances traveled by successful dispersers
	for(i in 1:numreefs){ # add up the immigrants from each patch reef
		j = abs(x)>(i*reefdist-reef) & abs(x)<(i*reefdist+reef)
		outside = outside + sum(j)
		dists = c(dists, x[j])
	}
	sum(abs(x)<reef)/(outside+sum(abs(x)<reef)) # 50% of settlers are local origin
	# reefdist = 10km and sd = 12km: 33% are local origin
	# reefdist = 15km and sd = 11km: 54% are local origin
	# reefdist = 15km and sd = 7km: 83% are local origin
	# reefdist = 10km and sd = 7km: 57% are local origin
	# reefdist = 10km and sd = 11km: 36% are local origin
	# reefdist = 15km and sd = 11km: 55% are local origin
	sd(x)
	sd(dists) # stdev of dispersal distances is 8km (or 12km)
	quantile(abs(dists), c(0.5, 0.75, 0.95))



#### 2 dimensional

	# Gaussian in a 2D continuous reef (retention = self-recruitment)
	len = 500000
	x = rnorm(len, sd=10)
	y = rnorm(len, sd=10)
	sum(abs(sqrt(x^2+y^2))<0.5)/len # 10km: only 0.12% retained on 1km diameter reef
	sum(abs(sqrt(x^2+y^2))<2)/len # 10km: only 1.9% retained on 4km diameter reef




######################################
###### Siegel et al 2003 MEPS ########


# Heuristic null model of dispersal

sigmaU = c(2.5,5,10,15, 20) # cm/s (stdev of alongshore current)
Tm = 11.5 # PLD in days
Tm = 11 # PLD in days # max from Thesher et al 1989
Tm = 7 # PLD in days # min from Thesher et al 1989
U = 0 # mean alongshore current

sigmaU = sigmaU * 3600*24/100/1000 # km/day

sigmaD = 2.238 * sigmaU * sqrt(Tm) # 33, 66, and 98 km
sigmaD # for Tm = 7, sigmaD = 25.6 km (sigmaU=5) or = 102 (sigmaU=20)
		# for Tm = 11, sigmaD = 32.1 km (sigmaU=5) or 128 (sigmaU=20)

Dd = 0.695 * Tm * U + 0.234 * Tm * sigmaU
Dd

Ddalt = sqrt(2/pi)*sigmaD # if xd = 0
Ddalt

mean(abs(rnorm(100000, mean=0, sd = sigmaD[1]))) # check Ddalt against a random number



##########################################
######## Doherty et al. 1995 #############
Can we get 4 migrants per generation across 1000km?

len = 1
migs = numeric(len)
for(i in 1:len){
	x = rnorm(n=1000, mean=0, sd=27) # n = # settlers per gen (=pop size for stable pop)
	migs[i] = sum(x>1000)	
}
max(migs)

len = 10000
x = rnorm(len, mean=0, sd = 7)
sum(abs(x)>=1000)/len
max(abs(x))
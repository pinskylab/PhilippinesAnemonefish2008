# Is sample He always greater than spp He?

x = 1:100 # sample He (very high!)
y = rep(1,100) # unsampled He (very low)


pop = x
pop = c(x,y)

len = 1000
u = round(runif(len, min=1, max=length(pop)))
v = round(runif(len, min=1, max=length(pop)))
match = 0
for(i in 1:len){
	match = match+(pop[u[i]] == pop[v[i]])
}
1-match/len
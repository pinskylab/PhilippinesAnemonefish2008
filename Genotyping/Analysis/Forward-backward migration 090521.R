# exploring forward and backward migration matrices

# Pop sizes
n = c(1,1,1) # even
n = c(10,1,1) # not even
	
# Forward migration matrix: f[i,j] is probability of going from i to j
# Rows must sum to 1
f = matrix(nrow = 3, data = c(
	1/3,1/3,1/3,
	1/3,1/3,1/3,
	1/3,1/3,1/3))

f = matrix(nrow = 3, data = c(
	0.1,0.45,0.45,
	0.45,0.1,0.45,
	0.45,0.45,0.1))

	
N = length(f)	
	
# Backward migration matrix: b[i,j] is fraction of i originating from j
# Rows sum to 1

b = matrix(nrow = nrow(f), data = rep(NA, N))

for(i in 1:nrow(b)){
	denom = sum(n*f[,i])
	b[i,] = n*f[,i]/denom
}
b
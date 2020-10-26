## Plot Gaussian dispersel kernels

# the hard way
r = rnorm(10000000)
plot(density(r))


# the exact way
sigma = 1
mu = 0

x = seq(-4,4,by=0.1)

y = 1/(sigma*sqrt(2*pi))*exp(-(x-mu)^2/(2*sigma^2))
#y2 = exp((-1/2)*x^2)/sqrt(2*pi)

plot(x,y,type="l", col=NA, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
polygon(x,y, col="grey", border=NA)
abline(v=1)


# the exact way:grey on blue
sigma = 1
mu = 0

x = seq(-4,4,by=0.1)

y = 1/(sigma*sqrt(2*pi))*exp(-(x-mu)^2/(2*sigma^2))
#y2 = exp((-1/2)*x^2)/sqrt(2*pi)

par(bg="dark blue")
plot(x,y,type="l", col=NA, xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
polygon(x,y, col="grey", border=NA)


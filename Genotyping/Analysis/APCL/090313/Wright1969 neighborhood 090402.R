# recreating figures on p. 304 in Wright 1969 Evolution and the Genetics of Populations Vol 2 The Theory of Gene Frequencies

x= seq(0, 10, by=0.1)
sigma=1
ynorm=(1/(sigma*sqrt(2*pi)))*exp(-x^2/(2*sigma^2))

a = 1/1
b= .5
b^(-2*a)*gamma(3*a)/gamma(a) # sigma2
y0 = (b^a)/(2*gamma(a+1)) 
ylepto = y0*exp(-b*x^(1/a))


plot(x, ynorm, type="l", ylim=c(0,.45))
lines(x,ylepto, col="red")



a = seq(1e-40,2, by=0.01)
kurt = gamma(a)*gamma(5*a)/(gamma(3*a)^2) -3 # kurtosis for the lepto curve
mult = 2^(a+1)*gamma(a+1)*(gamma(a)/gamma(3*a))^0.5 # multiplier for sigma*d to get N

plot(kurt, mult, type="l")
abline(v=0, lty=3)

hist(mult)
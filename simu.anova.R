# simulation of ANOVA main effects, 1-way
#

# \brief	compute moments from the data matrix
compute.moments	<- function(x) {
	x	= as.matrix(x)
	n	= nrow(x) # number of samples
	
	w 	= 1/n
	
	#m 	= apply(x, 2, mean)	
	m 	= w*colSums(x) # same as above, but much faster
	
	# v = apply(x, 2, var)
  	v 	= (n/(n-1)) * ( colSums(w*x^2)-colSums(w*x)^2 ) # same as above, but much faster

  	# set small values of variance exactly to zero
  	v[v < .Machine$double.eps] = 0
	
  	return( list(mean=m, var=v) )
}

generate_data	<- function(l, tau) {
	x	= rnorm(N, mean = nu, sd = tau)
	return(x)
}

simu.anova 		<- function(l) {
	# generate samples
	X 	= sapply(1:K, generate_data, tau)

	# calculate the statistical estimator
	moments = compute.moments(X)
	mu_		= moments$mean
	var_	= moments$var
	
	grand_mean	= sum(X) / (K*N)
	SST			= N*sum( (mu_ - grand_mean)^2 )
	MST			= SST / (K - 1)
	SSE			= (N - 1)*sum(var_)
	MSE			= SSE / (N*K - K)
	
	# F-statistic
	F 		= MST / MSE
	return (F)
}

# number of simulation
num		= 2000
# number of samples
N		= 3
# number of groups
K		= 5

require("geoR")
# generating
# the variance follows a scale-inverse-chi-square distr. scale-inv(d_0, s_0^2) with s_0^2 = 4 and d_0 = 4 (according Strimmer's)
d_0		= 4
s_0		= 2

tau		= sqrt( rinvchisq(1, df = d_0, scale = s_0) )

# the mean
nu		= runif(1, min = -10, max = 10)

x 		= sapply(1:num, simu.anova)
xx 		= sort(as.vector(x))

z 		= (1:num)/num
plot(xx, z, type = "l")
lines(xx, pf( xx, df1 = (K - 1), df2 = (N*K - K) ), col=3, lwd = 2)
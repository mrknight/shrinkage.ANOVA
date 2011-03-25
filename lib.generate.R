#################################
# \brief	generating data class for simulation of shrinkage ANOVA
#
#################################

# pre-defined variables

num		= 500 # number of simulation
P		= 2000 # number of genes
N		= 3 # sample size
K_		= c(3, 5, 7) # number of groups

# dir
dir.data = "/Users/knight/dev/data/shrinkage/"
dir.output = "/Users/knight/dev/output/shrinkage/"

require("geoR")
# generating
# the variance follows a scale-inverse-chi-square distr. scale-inv(d_0, s_0^2) with s_0^2 = 4 and d_0 = 4 (according to Strimmer's)
d_0		= 4
s_0		= 2
tau		= sqrt( rinvchisq(1, df = d_0, scale = s_0) )

# the mean
nu		= runif(1, min = -10, max = 10)
										
# \brief	generate data for the I. simulation (option I)
#			same mean, same variance
#			dumb variable l, p
generate_data_I	<- function(l, tau, K) {
	x	= rnorm(N*K, mean = nu, sd = tau)
	return(x)
}

# TODO
# \brief	generate data for the II. simulation (option II)
#			same mean, same variance
#			with p% outliers, variance of outliers = variance^2
#			dumb variable l
generate_data_II <- function(l, p) {
	if (runif(1, 0, 1) < p) {
		x	= rnorm(N, mean = nu, sd = tau^2)
	} else {
		x	= rnorm(N, mean = nu, sd = tau)
	}
	return(x)
}

# TODO:
simu <- function() {
	# generate the evaluated data, option I
	xd	= sapply(1:num, compute_statistic, data.info=optionI)

	# save evaluated data to file
	save(xd, file = "A.simuI.rdata")

	# generate the evaluated data, option III
	xd	= sapply(1:num, compute_statistic, data.info=optionII)

	# save evaluated data to file
	save(xd, file = "A.simuII.rdata")

}

K	= K_[2]
# generate samples
X 	= sapply(1:P, generate_data_I, tau, K)
return (X)


#XD = sapply(1:1, test)
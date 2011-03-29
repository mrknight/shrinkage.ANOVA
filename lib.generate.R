#################################
# \brief	generating data class for simulation of shrinkage ANOVA
#
#################################

# pre-defined variables

num		= 500 # number of simulation
P		= 5000 # number of genes
N		= 3 # sample size
K_		= c(2, 3, 5, 7) # number of groups
p_		= c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5) # vector of % outliers

# source dir
dir.data 	= "/Users/knight/dev/data/shrinkage.ANOVA"
dir.output 	= "/Users/knight/dev/output/shrinkage.ANOVA"
dir.source	= "/Users/knight/dev/src/shrinkage.ANOVA"

# loading my library
source(paste(dir.source, "/lib.evaluate.R"		, sep =""))
source(paste(dir.source, "/lib.shrinkvar.R"		, sep =""))
source(paste(dir.source, "/lib.output.R"		, sep =""))


#require("geoR")
# generating
# the variance follows a scale-inverse-chi-square distr. scale-inv(d_0, s_0^2) with s_0^2 = 4 and d_0 = 4 (according to Strimmer's)
#d_0		= 4; s_0		= 2
#tau		= sqrt( rinvchisq(1, df = d_0, scale = s_0) )
# the mean
#nu		= runif(1, min = -10, max = 10)

tau 		= 10
nu			= 0

# \brief	generate data from a normal dist with a given p
#			p was given as the number of outliers
#			same mean, same variance
#			dumb variable l
generate_data	<- function(l, p) {
	if (runif(1, 0, 1) < p) {
		x	= rnorm(N, mean = nu, sd = tau^2)
	} else {
		x	= rnorm(N, mean = nu, sd = tau)
	}
	return(x)
}
										
# \brief	generate data for the I. simulation (option I) without outliers
#			same mean, same variance
#			dumb variable l
generate_data_I	<- function(l, K) {
#	x	= rnorm(N*K, mean = nu, sd = tau)
	x	= sapply(1:K, generate_data, p = 0)
	return(x)
}

# \brief	generate data for the II. simulation (option II) with outlier
#			same mean, same variance
#			with p% outliers, variance of outliers = variance^2
#			dumb variable l
generate_data_II <- function(l, K, p) {
	x	= sapply(1:K, generate_data, p)
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

#optionI		= list(	simuart1 = "generate_data_I", simuart2 = "generate_data_I",
#					sigma1 = sigma1, sigma2 = sigma1, p = NA)
#optionII 	= list(	simuart1 = "generate_data_II", simuart2 = "generate_data_II",
#					sigma1 = sigma1, sigma2 = sigma1, p = 0.01)

# generate the evaluated data, option I
#xd	= sapply(1:num, compute_statistic, data.info=optionI)

# 
# \brief	testing the compute of lambda
#

num		= 100
# number of genes
P		= 2000
# number of samples
N		= 3

# \brief	compute moments from the data matrix
compute_moments <- function(x) {
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

# \brief	generate data for testing the lambda
#			normal dist. with mean = 0 and selectable variance (default var = 1)
generate_data	<- function(l, var = 1) {
	x	= rnorm(N, mean = 0, sd = sqrt(var))
	return(x)
}

# \brief	compute the shrinkage estimates with the hang over lambda, the targets 
#			and the observed data associated estimates
calc_shrinkage_var <- function(lambda, v_k, v_target) {
	v_k_ 	= lambda * v_target + (1 - lambda) * v_k
}

# \brief	calc the estimated pooling parameter lambda
#			
# \param	x as data
#			x_k_ as mean of each gene
#			v_k as unbiased variance of each gene
#			v_target as the target of the shrinkage method
#		
calc_lambda <- function(x, x_k_, v_k, v_target) {
	# number of genes
	k 		= ncol(x) # = P
	# samples size
	n 		= nrow(x) # = N
	# sample variance of each gene
	w_k_ 	= v_k * (n-1) / n	
	w_ik	= sweep(x, 2, x_k_, "-")^2

	# empirical variance of v_k Var^(v_k)
	var_vk 	= ( n / (n - 1)^3 ) * colSums( sweep(w_ik, 2, w_k_, "-")^2)
	numerator	= sum(var_vk)
	denominator = sum( (v_k - v_median)^2 )

	lambda 	= numerator/denominator
	lambda 	= min(1, lambda)
#	v_k_ 	= lambda * v_median + (1 - lambda) * v_k
	return(lambda)
}

xd		= sapply(1:N, generate_data)
#save(xd, file = "xd.rdata")

#load("xd.rdata")
x_m	 	= compute_moments(xd)
x_k_	= x_m$mean
v_k		= x_m$var
	
lambda	= calc_lambda(x, x_k_, v_k, v_median)

#############################################################################
# \brief	shrinkage variance class
#  			calculate the shrinkage variance estimator towards the median
#			contains all of methods
# \param	
#
# \return	
#############################################################################

# \brief	compute moments from the data matrix
compute.moments <- function(x) {
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

# \brief	compute the estimated pooling parameter lambda
#			
# \param	x 		: data
#			x_k_	: mean of each gene
#			v_k 	: unbiased variance of each gene
#			v_target: the target of the shrinkage method
#		
compute.lambda <- function(x, x_k_, v_k, v_target) {
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
	denominator = sum( (v_k - v_target)^2 )

	lambda 	= numerator/denominator
	lambda 	= min(1, lambda)
	
#	v_k_ 	= lambda * v_target + (1 - lambda) * v_k
	return(lambda)
}

# \brief	compute the shrinkage variance 
#			and the observed data associated estimates
compute.shrink.var	<- function(lambda, v_target, v_k) {
	sv 		= lambda * v_target + (1 - lambda) * v_k
	return (sv)
}

# \TODO		bugs?
#############################################################################
# \brief	compute F statistic as function of lambda
# \param
#
# \return	mean and variance of the function F
#############################################################################
compute.ab <- function(lambda, mean.all, var.all, target.all, K) {
	# shrinkage variances definition
	sv.all		= matrix(NA, nrow = K, ncol = P) 
	
	# compute all the shrinkage variances with a given lambda	
	for (i in 1:K) {
		sv.all[i,] 		= compute.shrink.var(lambda, target.all[i], var.all[i,])
	}
	
	F		= compute.F.score(mean.all, sv.all, K)
	
	a	= sum(F)
	b	= sum(F^2)
	return (list(a=a/P, b=b/(P-1)))
}

# \TODO		bugs?
#############################################################################
# \brief	compute F statistic as function of lambda for case 3, 5
# \param
#
# \return	mean and variance of the function F
#############################################################################
compute.ab2 <- function(lambda, mean.all, var.all, target, K) {
	
	# compute all the shrinkage variances with a given lambda	
	sv		= compute.shrink.var(lambda, target, var.all)
	# replicate the sv for K groups to compute the F statistics
	U.sv	= matrix(rep(sv, K), ncol = P, nrow = K, byrow = T)
	
	F		= compute.F.score(mean.all, U.sv, K)
	
	a	= sum(F)
	b	= sum(F^2)
	return (list(a=a/P, b=b/(P-1)))
}

# \TODO		not finish
# \brief	compute the score statistic of shrinkage variance from the data
# \param	X 	: data matrix
#			cl	: class label vector
#			var.equal : data with equal variance or not
# \notes	using the Strimmer's variant for choosing targets
stat.shrinkt <- function(X, cl, var.equal=TRUE) {
	if (var.equal) {
		
	}
#	(x, y) = data.process(X, L)
	xm = compute_moments(x)
	ym = compute_moments(y)
	
	xt = median(xm$var)
	yt = median(ym$var)
	lambda_x	= compute_lambda(x, xm$mean, xt)
	lambda_y	= compute_lambda(y, ym$mean, yt)
	
	v_k_ 	= calc_shrinkage_var(lambda_x, xm$var, xt)
	w_k_ 	= calc_shrinkage_var(lambda_y, ym$var, yt)
	
	#score   = (x_k - y_k) / sqrt(v_k_/L + w_k_/L)
}
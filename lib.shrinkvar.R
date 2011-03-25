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
#	v_k_ 	= lambda * v_median + (1 - lambda) * v_k
	return(lambda)
}

# \brief	compute the shrinkage variance 
#			and the observed data associated estimates
shrink.var	<- function(x) {
	
	sv 		= lambda * v_target + (1 - lambda) * v_k	
}

# \brief	generate data for testing the lambda
#			normal dist. with mean = 0 and selectable variance (default var = 1)
#			l as dumb variable
generate_data	<- function(l, var = 1) {
	x	= rnorm(N, mean = 0, sd = sqrt(var))
	return(x)
}


# \brief	compute t statistic as function of lambda
# \param
#
# \return
compute_t <- function(lambda, param) {
	if (is.na(param$w_k)[1]) {
		v_k_ 	= calc_shrinkage_var(lambda, param$v_k, param$v_target)
		t_k		= (param$x_k - param$y_k) / sqrt(2*v_k_/L)		
	}
	else {
		v_k_ 	= calc_shrinkage_var(lambda, param$v_k, param$v_target)
		w_k_ 	= calc_shrinkage_var(lambda, param$w_k, param$w_target)	
		t_k		= (param$x_k - param$y_k) / sqrt(v_k_/L + w_k_/L)
	}
	a	= sum(t_k)
	b	= sum(t_k^2)
	return (list(a=a/N, b=b/(N-1)))
}

# \brief	compute all the variations of scores and the corresponded t statistic
compute_statistic <- function(l, data.info) {
	# without centering
	center 	= FALSE
	x		= sapply(1:N, data.info$simuart1, data.info$p)
	y		= sapply(1:N, data.info$simuart2, data.info$p)
	
	sigma1 	= data.info$sigma1
	sigma2 	= data.info$sigma2	
	sigma  	= (sigma1 + sigma2) / 2
	
	# moments of x
	x_m	 	= compute_moments(x)
	x_k		= x_m$mean
	v_k		= x_m$var
	
	# moments of y
	y_m	 	= compute_moments(y)
	y_k		= y_m$mean
	w_k		= y_m$var

	# 
	lambda	= seq(0.01, .99, length=99)
	# 1. variant (the same with Strimmer's if the variances in each group are different)
	v_median 	= median(v_k)
	w_median	= median(w_k)
	lambda_v	= compute_lambda(x, x_k, v_median)
	lambda_w	= compute_lambda(y, y_k, w_median)
	lambda1		= c(lambda_v, lambda_w)
	
	v_k_ 	= calc_shrinkage_var(lambda_v, v_k, v_median)
	w_k_ 	= calc_shrinkage_var(lambda_w, w_k, w_median)
	u_k_ 	= (v_k_ + w_k_) / 2
		
	score_1	 = sum((v_k_ - sigma1^2)^2) + sum((w_k_ - sigma2^2)^2)
	scoret_1 = sum((u_k_ - sigma^2)^2)
	t_k1 	 = (x_k - y_k) / sqrt(v_k_/L + w_k_/L)

	# function of t as lambda
	param = list(x_k=x_k, y_k=y_k, v_k=v_k, w_k=w_k, v_target=v_median, w_target=w_median)	
	tt = sapply(lambda, compute_t, param)
	a1 = as.vector(tt[1,],"numeric")
	b1 = as.vector(tt[2,],"numeric")	
	
	# 2. variant
	x_median	= median(c(v_k, w_k))
	lambda_v	= compute_lambda(x, x_k, x_median)
	lambda_w	= compute_lambda(y, y_k, x_median)
	lambda2		= c(lambda_v, lambda_w)

	v_k_ 	= calc_shrinkage_var(lambda_v, v_k, x_median)
	w_k_ 	= calc_shrinkage_var(lambda_w, w_k, x_median)
	u_k_ 	= (v_k_ + w_k_) / 2

	score_2	 = sum((v_k_ - sigma1^2)^2) + sum((w_k_ - sigma2^2)^2)
	scoret_2 = sum((u_k_ - sigma^2)^2)
	t_k2 	 = (x_k - y_k) / sqrt(v_k_/L + w_k_/L)
	
	param = list(x_k=x_k, y_k=y_k, v_k=v_k, w_k=w_k, v_target=x_median, w_target=x_median)	
	tt = sapply(lambda, compute_t, param)
	a2 = as.vector(tt[1,],"numeric")
	b2 = as.vector(tt[2,],"numeric")
	
	# 3. variant
	u_k			= (v_k + w_k) / 2
	u_median	= median(u_k)
	
	# combine x and y
	xy			= matrix(0, ncol=N, nrow=2*L)
	xy[1:L, 1:N] = x[1:L, 1:N]
	xy[(L+1):(2*L), 1:N] = y[1:L, 1:N]
	xy_m	= compute_moments(xy)
	
	# TODO : general check ???
	lambda_u	= compute_lambda(xy, xy_m$mean, u_median)
	lambda3		= c(lambda_u)

	u_k_ 	= calc_shrinkage_var(lambda_u, u_k, u_median)
	
	score_3	 = NA
	scoret_3 = sum((u_k_ - sigma^2)^2)
	t_k3 	 = (x_k - y_k) / sqrt(2*u_k_/L)
	
	param = list(x_k=x_k, y_k=y_k, v_k=u_k, w_k=NA, v_target=u_median, w_target=NA)
	tt = sapply(lambda, compute_t, param)
	a3 = as.vector(tt[1,],"numeric")
	b3 = as.vector(tt[2,],"numeric")
	
	# 4. variant
	u_k			= (v_k + w_k) / 2
	u_median	= median(u_k)
	
	lambda_v	= compute_lambda(x, x_k, u_median)
	lambda_w	= compute_lambda(y, y_k, u_median)
	lambda4		= c(lambda_v, lambda_w)
	
	v_k_ 	= calc_shrinkage_var(lambda_v, v_k, u_median)
	w_k_ 	= calc_shrinkage_var(lambda_w, w_k, u_median)
	u_k_ 	= (v_k_ + w_k_) / 2

	score_4	 = sum((v_k_ - sigma1^2)^2) + sum((w_k_ - sigma2^2)^2)
	scoret_4 = sum((u_k_ - sigma^2)^2)
	t_k4 	 = (x_k - y_k) / sqrt(v_k_/L + w_k_/L)
	
	#function of t as lambda
	param = list(x_k=x_k, y_k=y_k, v_k=v_k, w_k=w_k, v_target=u_median, w_target=u_median)
	tt = sapply(lambda, compute_t, param)
	a4 = as.vector(tt[1,],"numeric")
	b4 = as.vector(tt[2,],"numeric")
		
	# 5. variant (from Strimmer, without centering), only for data with the same variance

	tgt		 = median(xy_m$var)
	lambda_u = compute_lambda(xy, xy_m$mean, tgt)
	u_k_ 	 = calc_shrinkage_var(lambda_u, xy_m$var, tgt)
	lambda5  = c(lambda_u)

	score_5	 = NA
	scoret_5 = sum((u_k_ - sigma^2)^2)
	t_k5 	 = (x_k - y_k) / sqrt(2*u_k_/L)
	
	#function of t as lambda
	param = list(x_k=x_k, y_k=y_k, v_k=xy_m$var, w_k=NA, v_target=tgt, w_target=NA)
	tt = sapply(lambda, compute_t, param)
	a5 = as.vector(tt[1,],"numeric")
	b5 = as.vector(tt[2,],"numeric")
	
	return( list(score_1 = score_1, scoret_1 = scoret_1, t1 = t_k1, a1 = a1, b1 = b1,
 				 score_2 = score_2, scoret_2 = scoret_2, t2 = t_k2, a2 = a2, b2 = b2,
				 score_3 = score_3, scoret_3 = scoret_3, t3 = t_k3, a3 = a3, b3 = b3,
				 score_4 = score_4, scoret_4 = scoret_4, t4 = t_k4, a4 = a4, b4 = b4,
				 score_5 = score_5, scoret_5 = scoret_5, t5 = t_k5, a5 = a5, b5 = b5,
				 lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4, lambda5 = lambda5) )
}

# \brief	
data.process <- function(X, L) {
	
}

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
# 
# \brief	testing the compute of lambda
#

num		= 200
# number of genes
N		= 2000
# number of samples
L		= 3


# \brief	compute the shrinkage estimates with the hang over targets 
#			and the observed data associated estimates
calc_shrinkage_var <- function(lambda, v_k, v_target) {
	v_k_ 	= lambda * v_target + (1 - lambda) * v_k
}

# \brief	generate data for testing the lambda
#			normal dist. with variance = position of genes
generate_data	<- function(var) {
	x	= rnorm(L, mean = 0, sd = sqrt(var))
	return(x)
}


generate_data_1	<- function(sd) {
	x	= rnorm(L, mean = 0, sd = 1)
	return(x)
}

# \brief	compute moments
compute_moments <- function(x) {
	x	= as.matrix(x)
	n	= nrow(x)
	
	w 	= 1/n
	
	#m 	= apply(x, 2, mean)	
	m 	= colSums(w*x) # same as above, but much faster
	
	# v = apply(x, 2, var)
  	v 	= (n/(n-1)) * (colSums(w*x^2)-colSums(w*x)^2) # same as above, but much faster

  	# set small values of variance exactly to zero
  	v[v < .Machine$double.eps] = 0
	
  	return( list(mean=m, var=v) )
}

# \brief	get lambda from Strimmer's library
get_lambda <- function(x, mean, target) {
	# centre the data
	xc		= sweep(x, 2, mean)
#	xc 		= x
	n		= nrow(x)
	w 		= rep(1/n, n)
	
	z 		= pvt.get.lambda(xc, -1, w, verbose=TRUE, type="variance", target)
	
	return(z$lambda.var)
}

# \brief
test_lambda <- function(l) {
	x		= sapply(1:N, generate_data)
	# moments of x
	x_m	 	= compute_moments(x)
	x_k		= x_m$mean
	v_k		= x_m$var
	
	v_median 	= median(v_k)
	lambda_v	= get_lambda(x, x_k, v_median)
#	lambda_v	= calc_lambda(x, x_k, v_k, v_median)

	
#	xc 			= sweep(x, 2, x_k, "-")
#	lambda_v 	= calc_lambda(xc, compute_moments(xc)$mean, compute_moments(xc)$var, median(compute_moments(xc)$var))
	
	return(lambda_v)
}

# \brief
calc_shrinkage_t <- function(l) {
	x		= sapply(1:N, generate_data)
	y		= sapply(1:N, generate_data)
	# moments of x
	x_m	 	= compute_moments(x)
	x_k		= x_m$mean
	v_k		= x_m$var
	
	# moments of y
	y_m	 	= compute_moments(y)
	y_k		= y_m$mean
	w_k		= y_m$var
	
	#u_k		= (v_k + w_k) / 2
	v_median 	= median(v_k)
	w_median	= median(w_k)
	#u_median	= median(u_k)
	#x_median	= median(c(v_k, w_k))
	lambda_v	= get_lambda(x, x_k, v_median)
	#xc 			= sweep(x, 2, x_k, "-")
	#lambda_v1 	= calc_lambda(xc, x_k, compute_moments(xc)$var, median(compute_moments(xc)$var))
	
	#lambda_w = calc_lambda(y, y_k, w_k, w_median)	
	lambda_w	= get_lambda(y, y_k, w_median)
	
	v_k_ 	= calc_shrinkage_var(lambda_v, v_k, v_median)
	w_k_ 	= calc_shrinkage_var(lambda_w, w_k, w_median)
	
	t_k = (x_k - y_k) / sqrt(v_k_/L + w_k_/L)
	return(t_k)
}

# \brief	get shrinkage t statistic from Strimmer
get_shrinkage_t <- function(l) {
	x		= sapply(1:N, generate_data_1)
	y		= sapply(1:N, generate_data_1)
	X		= matrix(0, ncol=N, nrow=2*L)
	X[1:L, 1:N] = x[1:L, 1:N]
	X[(L+1):(2*L), 1:N] = y[1:L, 1:N]
	label	= c(rep(1, L), rep(2, L))
	
	t 		= shrinkt.stat(X, label, verbose = FALSE)
	return (t)
}

# \brief
compare_shrinkage_t <- function(l) {
	x		= sapply(1:N, generate_data_1)
	y		= sapply(1:N, generate_data_1)
	# moments of x
	x_m	 	= compute_moments(x)
	x_k		= x_m$mean
	v_k		= x_m$var
	
	# moments of y
	y_m	 	= compute_moments(y)
	y_k		= y_m$mean
	w_k		= y_m$var
	
	v_median 	= median(v_k)
	w_median	= median(w_k)
	#lambda_v = calc_lambda(x, x_k, v_k, v_median)
	lambda_v	= get_lambda(x, x_k, v_median)
	#lambda_w = calc_lambda(y, y_k, w_k, w_median)	
	lambda_w	= get_lambda(y, y_k, w_median)
	
	v_k_ 	= calc_shrinkage_var(lambda_v, v_k, v_median)
	w_k_ 	= calc_shrinkage_var(lambda_w, w_k, w_median)
	
	diff1	= (x_k - y_k)
	sd1		= sqrt(v_k_/N + w_k_/N) 
	t_k =  diff1/sd1 
	
	X		= matrix(0, ncol=N, nrow=2*L)
	X[1:L, 1:N] = x[1:L, 1:N]
	X[(L+1):(2*L), 1:N] = y[1:L, 1:N]
	label	= c(rep(1, L), rep(2, L))
	
	tmp = centroids(X, label, var.pooled=TRUE, var.groups=FALSE, shrink=TRUE, verbose=TRUE)

    # differences between the two groups
    diff = tmp$means[,1]-tmp$means[,2]

    # standard error of diff
    n1 = tmp$samples[1]
    n2 = tmp$samples[2]
    v =  tmp$var.pooled   
    sd = sqrt( (1/n1 + 1/n2)*v )

	#t 		= shrinkt.stat(X, label, verbose = FALSE)
	
	xc		= sweep(x, 2, x_k, "-")
	yc		= sweep(y, 2, y_k, "-")
	XC		= matrix(0, ncol=N, nrow=2*L)
	XC[1:L, 1:N] = xc[1:L, 1:N]
	XC[(L+1):(2*L), 1:N] = yc[1:L, 1:N]
	
	
#	print(t)
#	print(t_k)
#	return(t_k)
}

#yy	= (1:num)/num
#df 	= L + L
xd	= sapply(1:num, test_lambda)
hist(xd)
#xx  = sort(xd[1,])
#plot(xx, yy, type="l");
#lines(xx, pt(xx, df), col=4, lwd = 1)

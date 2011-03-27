#################################
#
# \brief	evaluate the given data
#
#################################

# \brief	convert 1 column of X (data of 1 gene) to NxK matrix 
convert.col2matrix 	<- function(x.col, K) {
	return (matrix(x.col, nrow = N, ncol = K))
}

# compute the mean for each genes, each groups
compute.mean.all <- function(pos, X, K) {
	x.col	= X[,pos]
	X.col 	= convert.col2matrix(x.col, K)
	moments = compute.moments(X.col)
	return(moments$mean)
}

# compute the var for each genes, each groups	
compute.var.all <- function(pos, X, K) {
	x.col	= X[,pos]
	X.col 	= convert.col2matrix(x.col, K)
	moments = compute.moments(X.col)
	return(moments$var)
}

# \brief	evalutate the F statistic
eval.statistic <- function(X, K, option) {
	P	= ncol(X) # number of genes
	N	= nrow(X) / K # number of samples
	
	# \brief	compute ordinary F statistic
	compute.F <- function(X = X.col) {
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
	compute.F.all <- function(pos, X) {
		x.col	= X[,pos]
		X.col 	= convert.col2matrix(x.col)
		F		= compute.F(X.col)
		return(F)
	}
	
	compute.shrinkage.f <- function() {
		
	}
	
	f.all = sapply(1:P, compute.F.all, X)
	return (f.all)
}

# \brief	evalutate the F statistic
compute.F.score <- function(mean.all, var.all, K) {
	grandmean.all	= colSums(mean.all) / K
	
	SST			= N*colSums( t(t(mean.all) - grandmean.all)^2 )
	MST			= SST / (K - 1)
	SSE			= (N - 1)*colSums(var.all)
	MSE			= SSE / (N*K - K)

	# F-statistic
	F 		= MST / MSE	
	return (F)
}

# \brief	compute all the variations of scores and the corresponded F statistic
#			generate and evaluate data on the fly
eval.F.statistic <- function(l, K, p) {
	# header definition
	lambda.all	= rep(0, K)
	sv.all		= matrix(NA, nrow = K, ncol = P) 
	
	# without centering
	center 	= FALSE
	
	sigma 		= tau
	
	# generate the simulated data
	X 	= sapply(1:P, generate_data_II, K, p)

	mean.all 	= sapply(1:P, compute.mean.all, X, K)
	var.all		= sapply(1:P, compute.var.all, X, K)
	grandmean.all	= colSums(mean.all) / K

	# 1.case (the same with Strimmer's if the variances in each group are different)
	target.all	= apply(var.all, 1, median)
	
	for (i in 1:K) {
		lambda.all[i]	= compute.lambda(X[(1+N*(i-1)):(N*i),], mean.all[i,], var.all[i,], target.all[i])
		sv.all[i,] 		= compute.shrink.var(lambda.all[i], target.all[i], var.all[i,])
	}
	
	lambda1		= lambda.all
	
	u.all 		= colMeans(sv.all)
		
	score_1	 	= sum((sv.all - sigma^2)^2)
	scoret_1 	= sum((u.all - sigma^2)^2)
	F.all_1 	= compute.F.score(mean.all, sv.all, K)
	
	# 2.case
	target		= median(var.all)
	for (i in 1:K) {
		lambda.all[i]	= compute.lambda(X[(1+N*(i-1)):(N*i),], mean.all[i,], var.all[i,], target)
		sv.all[i,] 		= compute.shrink.var(lambda.all[i], target, var.all[i,])
	}
	
	lambda2		= lambda.all

	u.all 		= colMeans(sv.all)
		
	score_2	 	= sum((sv.all - sigma^2)^2)
	scoret_2 	= sum((u.all - sigma^2)^2)
	F.all_2 	= compute.F.score(mean.all, sv.all, K)
		
	# 3.case
	u.all		= colMeans(var.all)
	u.median	= median(u.all)
	
	# combined moments of all groups
	# treat all groups like one group
	X.m			= compute.moments(X)

	lambda3		= compute.lambda(X, X.m$mean, X.m$var, u.median)

	# compute the shrinkage variance
	u.sv		= compute.shrink.var(lambda3, u.median, X.m$var)
	
	score_3	 	= NA
	scoret_3 	= sum((u.sv - sigma^2)^2)
	U.sv		= matrix(rep(u.sv, K), ncol = P, nrow = K, byrow = T)
	F.all_3  	= compute.F.score(mean.all, U.sv, K)
	
	#param = list(x_k=x_k, y_k=y_k, v_k=u_k, w_k=NA, v_target=u_median, w_target=NA)
	#tt = sapply(lambda, compute_t, param)
	#a3 = as.vector(tt[1,],"numeric")
	#b3 = as.vector(tt[2,],"numeric")
	
	# 4.case
	u.all		= colMeans(var.all)
	u.median	= median(u.all)
	
	for (i in 1:K) {
		lambda.all[i]	= compute.lambda(X[(1+N*(i-1)):(N*i),], mean.all[i,], var.all[i,], u.median)
		sv.all[i,] 		= compute.shrink.var(lambda.all[i], u.median, var.all[i,])
	}
	
	lambda4		= lambda.all
	u.all 		= colMeans(sv.all)
		
	score_4	 	= sum((sv.all - sigma^2)^2)
	scoret_4 	= sum((u.all - sigma^2)^2)
	F.all_4 	= compute.F.score(mean.all, sv.all, K)
	
	# 5.case (aslike from Strimmer,  but without centering), only for data with the same variance

	z.median	= median(X.m$var)
	lambda5 	= compute.lambda(X, X.m$mean, X.m$var, z.median)
	
	# compute the shrinkage variance
	z.sv		= compute.shrink.var(lambda5, z.median, X.m$var)
	
	score_5	 	= NA
	scoret_5 	= sum((z.sv - sigma^2)^2)
	Z.sv		= matrix(rep(z.sv, K), ncol = P, nrow = K, byrow = T)
	F.all_5  	= compute.F.score(mean.all, Z.sv, K)

	# dumb variable 
	# TODO: add a and b
	a1 = a2 = a3 = a4 = a5 = 0
	b1 = b2 = b3 = b4 = b5 = 0
	
	return( list(score_1 = score_1, scoret_1 = scoret_1, F_1 = F.all_1, a1 = a1, b1 = b1,
 				 score_2 = score_2, scoret_2 = scoret_2, F_2 = F.all_2, a2 = a2, b2 = b2,
				 score_3 = score_3, scoret_3 = scoret_3, F_3 = F.all_3, a3 = a3, b3 = b3,
				 score_4 = score_4, scoret_4 = scoret_4, F_4 = F.all_4, a4 = a4, b4 = b4,
				 score_5 = score_5, scoret_5 = scoret_5, F_5 = F.all_5, a5 = a5, b5 = b5,
				 lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4, lambda5 = lambda5) )
}

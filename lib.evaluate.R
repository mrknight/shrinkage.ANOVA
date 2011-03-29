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

# \TODO		obsolete?
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

# \brief	compute the F statistic with given means and variances
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
	# variables definition
	lambda.all	= rep(0, K)
	sv.all		= matrix(NA, nrow = K, ncol = P) 
	sigma 		= tau
	lambda.vec	= seq(0, 1, length=101)
	# without centering
	center 	= FALSE
	
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
	
	# F as function of lambda
	f 	= sapply(lambda.vec, compute.ab, mean.all, var.all, target.all, K)
	a1 	= as.vector(f[1,],"numeric")
	b1 	= as.vector(f[2,],"numeric")
	
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
	
	f 	= sapply(lambda.vec, compute.ab, mean.all, var.all, rep(target, K), K)
	a2 	= as.vector(f[1,],"numeric")
	b2 	= as.vector(f[2,],"numeric")
	
	# 3.case
	u.all		= colMeans(var.all)
	u.median	= median(u.all)
	
	# combined moments of all groups
	# treat all groups like one group
	X.m			= compute.moments(X)

	lambda3		= compute.lambda(X, X.m$mean, X.m$var, u.median)

	# compute the shrinkage variance
	u.sv		= compute.shrink.var(lambda3, u.median, u.all)
	
	score_3	 	= NA
	scoret_3 	= sum((u.sv - sigma^2)^2)
	# replicate the sv for K groups to compute the F statistics
	U.sv		= matrix(rep(u.sv, K), ncol = P, nrow = K, byrow = T)
	F.all_3  	= compute.F.score(mean.all, U.sv, K)
	
	f 	= sapply(lambda.vec, compute.ab2, mean.all, u.all, u.median, K)
	a3 	= as.vector(f[1,],"numeric")
	b3 	= as.vector(f[2,],"numeric")
	
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
	
	f 	= sapply(lambda.vec, compute.ab, mean.all, var.all, rep(u.median, K), K)
	a4 	= as.vector(f[1,],"numeric")
	b4 	= as.vector(f[2,],"numeric")
	
	# 5.case (like from Strimmer,  but without centering), only for data with the same variance (true with the null hypothesis of ANOVA)
	z.median	= median(X.m$var)
	lambda5 	= compute.lambda(X, X.m$mean, X.m$var, z.median)
	
	# compute the shrinkage variance
	z.sv		= compute.shrink.var(lambda5, z.median, X.m$var)
	
	score_5	 	= NA
	scoret_5 	= sum((z.sv - sigma^2)^2)
	Z.sv		= matrix(rep(z.sv, K), ncol = P, nrow = K, byrow = T)
	F.all_5  	= compute.F.score(mean.all, Z.sv, K)

	f 	= sapply(lambda.vec, compute.ab2, mean.all, X.m$var, z.median, K)
	a5 	= as.vector(f[1,],"numeric")
	b5 	= as.vector(f[2,],"numeric")
	
	return( list(score_1 = score_1, scoret_1 = scoret_1, F_1 = F.all_1, a1 = a1, b1 = b1,
 				 score_2 = score_2, scoret_2 = scoret_2, F_2 = F.all_2, a2 = a2, b2 = b2,
				 score_3 = score_3, scoret_3 = scoret_3, F_3 = F.all_3, a3 = a3, b3 = b3,
				 score_4 = score_4, scoret_4 = scoret_4, F_4 = F.all_4, a4 = a4, b4 = b4,
				 score_5 = score_5, scoret_5 = scoret_5, F_5 = F.all_5, a5 = a5, b5 = b5,
				 lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, lambda4 = lambda4, lambda5 = lambda5) )
}

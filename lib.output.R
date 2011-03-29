#################################
#
# \brief	library for output methods
#
#################################

find_minmax <- function(xd, data_pos) {
	c_max = c()
	c_min = c()
	for (i in 1:(length(data_pos))) {
		max = max(as.vector(xd[data_pos[i],],"numeric"))
		min = min(as.vector(xd[data_pos[i],],"numeric"))
		c_max = c(c_max, max)
		c_min = c(c_min, min)
	}
	x_max = max(c_max)
	x_min = min(c_min)
	return (list(x_max=x_max, x_min=x_min))
}

plot_score <- function(xd, data_pos, info) {
	x_max = find_minmax(xd, data_pos)$x_max
	x_min = find_minmax(xd, data_pos)$x_min
	plot(NULL, NULL, xlim=c(x_min,x_max), ylim=c(0,1), main = info, xlab = "score", ylab = "Density", log = "x")
	lg = c(); cl = c()
	for (i in 1:(length(data_pos))) {
		# calc which variant from the position of data in list
		variant 	= (data_pos[i] %/% 5) + 1
		xx 		= sort(as.vector(xd[data_pos[i],],"numeric"))
		num		= length(xx)
		z 		= (1:num)/num

		lines(xx, z, col=cols[variant], lwd = 2, xlog = TRUE)		
		lg = c(lg, paste(variant ,". variant"))
		cl = c(cl, cols[variant])
	}

	legend(x="bottomright",inset=0.01,legend=lg,col=cl,lwd=5)
}

plot_score_t <- function(xd, data_pos, info) {
	x_max = find_minmax(xd, data_pos)$x_max
	x_min = find_minmax(xd, data_pos)$x_min

	lg = c(); cl = c()

	variant 	= (data_pos[1] %/% 5) + 1

	xx 		= sort(as.vector(xd[data_pos[1],],"numeric"))
	num		= length(xx)
	z 		= (1:num)/num

	plot(xx, z, xlim=c(x_min,x_max), ylim=c(0,1), type ="l", main = info, xlab = "scoret", ylab = "Density", log = "x", col=cols[variant], lwd = 2, cex.main = 2)
	
	lg = c(lg, paste(variant ,". case"))
	cl = c(cl, cols[variant])
	for (i in 2:(length(data_pos))) {
		# calc which variant from the position of data in list
		variant 	= (data_pos[i] %/% 5) + 1

		xx 		= sort(as.vector(xd[data_pos[i],],"numeric"))
		num		= length(xx)
		z 		= (1:num)/num

		lines(xx, z, col=cols[variant], lwd = 2, xlog = TRUE)		
		lg = c(lg, paste(variant ,". case"))
		cl = c(cl, cols[variant])
	}
	
	legend(x="bottomright",inset=0.01,legend=lg,col=cl,lwd=2)
}

plot_all_score <- function(K, p) {
	load(paste(dir.data,	"/K",K,"_p",p,".rdata", sep=""))

	data_pos = c(1, 6, 16)
	plot_score(ed, data_pos, info = paste("K = ",K,"; p = ",p,sep=""))
}

plot_all_score_t <- function(K, p) {
	load(paste(dir.data,	"/K",K,"_p",p,".rdata", sep=""))

	data_pos = c(2, 7, 12, 17, 22)
	plot_score_t(ed, data_pos, info = paste("K = ",K,"; p = ",p,sep=""))
}

output.F.data <- function(data, K, p) {
	# \TODO		check ?
	plot_ab <- function(data_pos) {
		# compute the mean of all a and b from all simulation
		a	= apply(do.call("cbind", data[data_pos,]), 1, mean)
		b	= apply(do.call("cbind", data[(data_pos+1),]), 1, mean)
	  	a[a < .Machine$double.eps] = 0

		# compute the dof d2 from a (mean of F dist)
		d2 	= (mean(a)*2) / (mean(a) - 1)

		# compute the dof d1 from b and d2 (var of F dist)
		d1	= 2*(d2^2)*(d2-2) / ( mean(b)*(d2-2)^2*(d2-4) - 2*(d2^2) )

		return (list(d1=d1, d2=d2))
	}

	data_pos 	= c(4, 9, 14, 19, 24)
	dof.lam		= sapply(data_pos, plot_ab)

	#pdf("a.pdf")
	yy		= (1:(num*P))/(num*P)
	# defree of freedom, assumed
	dof.theo.d1 = (K - 1)
	dof.theo.d2 = (N*K - K)
	
	png(filename=paste(dir.output,"/plotf_K",K,"_p",p,".png",sep=""), width=1700, height=2400, pointsize=15)
	par(mfrow=c(5,2))
	
	plot.Fstat <- function(data_pos) {
		case = (data_pos %/% 5) + 1	
		# turn list of f-statistic to matrix
		xts	= do.call("cbind", data[data_pos,])
		xx  = sort(as.vector(xts)) # for all genes, or xts[1,] for just the 1.gene	
		plot(xx, yy, type="l", lwd = 1, main = paste("K = ",K,",p = ",p,", ",case,".case",sep=""), xlab = "F-statistic", ylab = "Density", col = 1, cex.main = 2, cex.axis = 2, cex.lab = 2);
		# mean and variance of all F
		mean.F	= mean(xx)
		var.F 	= var(xx)
		
		# compute the dof d2 from a (mean of F dist)
		dof.var.d2 	= 2*mean.F/(mean.F - 1)
		# temporal variable
		d2 	= dof.var.d2
		# compute the dof d1 from b and d2 (var of F dist)
		dof.var.d1	= 2*(d2^2)*(d2-2) / ( var.F*(d2-2)^2*(d2-4) - 2*(d2^2) )
			
		lines(xx, pf(xx, df1 = dof.var.d1, df2 = dof.var.d2), col=2, lwd=1)
		lines(xx, pf(xx, df1 = dof.lam[,case]$d1, df2 = dof.lam[,case]$d2), col=3, lwd=1)
		lines(xx, pf(xx, df1 = dof.theo.d1, df2 = dof.theo.d2), col=4, lwd=1)

		lg = c(	"Cdf of shrinkage-F",
				paste("Cdf of F with estimated dof from variance and mean, d1 = "	,round(dof.var.d1,digit=2),			",d2 =",round(dof.var.d2,digit=2)			,sep=""),
				paste("Cdf of F with estimated dof from lambda, d1 = "				,round(dof.lam[,case]$d1,digit=2),	",d2 =",round(dof.lam[,case]$d2,digit=2)	,sep=""), 
				paste("Cdf of F with theoretical dof, d1 = "						,round(dof.theo.d1,digit=2),		",d2 =",round(dof.theo.d2,digit=2)			,sep=""))
		cl = c(1,2,3,4)
		plot(xx, yy-pf(xx, df1 = dof.var.d1, df2 = dof.var.d2), type="l", col=2, lwd=1, main = paste("Magnification for K = ",K,",p = ",p,", ",case,".case",sep=""), 
			xlab = "F-statistic", ylab = "Density", ylim = c(-.05, .05), cex.main = 2, cex.axis = 2, cex.lab = 2)
		lines(xx, yy-pf(xx, df1 = dof.lam[,case]$d1, df2 = dof.lam[,case]$d2), col=3, lwd=1)
		lines(xx, yy-pf(xx, df1 = dof.theo.d1, df2 = dof.theo.d2), col=4, lwd=1)
		abline(h = 0, col=1, lwd=2)

		legend(x="bottomright",inset=0.001,legend=lg,col=cl,lwd=5, cex = 2)

	#	ks = ks.test(yy, pt(xx, df))$statistic
	#	ks2 = ks.test(yy, pt(xx, df2))$statistic
	#	ks_dof = ks.test(yy, pt(xx, dof[variant]))$statistic
	#	cat(ks," ",ks2," ",ks_dof,"\n")
	}

	data_pos = c(3, 8, 13, 18, 23)
	sapply(data_pos, plot.Fstat)

}

eval_data2 <- function(data, K, option) {

	find_minmax <- function(data_pos) {
		c_max = c()
		c_min = c()
		for (i in 1:(length(data_pos))) {
			max = max(as.vector(data[data_pos[i],],"numeric"))
			min = min(as.vector(data[data_pos[i],],"numeric"))
			c_max = c(c_max, max)
			c_min = c(c_min, min)
		}
		x_max = max(c_max)
		x_min = min(c_min)
		return (list(x_max=x_max, x_min=x_min))
	}

	cols = rainbow(3*2)
	
	png(filename=paste("hist_score",option,".png",sep=""), width=2000, height=1000, pointsize=12)
	#par(mfrow=c(2,2))
	# compare the score
	plot_score <- function(data_pos) {
		x_max = find_minmax(data_pos)$x_max
		x_min = find_minmax(data_pos)$x_min
		plot(NULL, NULL, xlim=c(x_min,x_max), ylim=c(0,1), main = paste("Option ",option,sep=""), xlab = "score", ylab = "Density")
		lg = c()
		for (i in 1:(length(data_pos))) {
			# calc which variant from the position of data in list
			variant 	= (data_pos[i] %/% 5) + 1
			h 			= hist(as.vector(data[data_pos[i],],"numeric"), breaks = 100, plot = FALSE )
			h$density 	= h$counts/sum(h$counts)
			h$counts    = cumsum(h$counts)
			h$density   = cumsum(h$density)

			lines(h$mid, h$density, col=cols[i], lwd = 2)
			lg = c(lg, paste(variant ,". variant"))
		}

		legend(x="bottomright",inset=0.01,legend=lg,col=cols,lwd=3)
	}
	data_pos = c(1, 6, 16)
	plot_score(data_pos)
	
	png(filename=paste("hist_scoret",option,".png",sep=""), width=2000, height=1000, pointsize=12)
	# compare the score~
	data_pos = c(2, 7, 12, 17, 22)
	plot_score(data_pos)
	
	png(filename=paste("plot_a_b",option,".png",sep=""), width=1500, height=2500, pointsize=12)
	par(mfrow=c(5,3))
	# compare the expected value and variance of t-statistic
	dof = rep(1, 5)	
	lambda	= seq(0.01, .99, length=99)
	
	plot_ab <- function(data_pos) {
		# calc which variant from the position of data in list
		variant = (data_pos %/% 5) + 1
		
		# compute the mean of all a and b from all simulation
		a	= apply(do.call("cbind", data[data_pos,]), 1, mean)
		b	= apply(do.call("cbind", data[(data_pos+1),]), 1, mean)
	  	a[a < .Machine$double.eps] = 0
		# compute the defree of freedom from b
		dof[variant] = (mean(b)*2) / (mean(b) - 1)

		plot(lambda, a, type="l", lwd = 1, main = paste("Option ", option ,", ", variant ,". variant",sep=""), ylab = "Expected Value", xlab = "lambda")
		lines(lambda, rep(0, 99), col=4, lwd = 1)
		plot(lambda, b, type="l", lwd = 1, main = paste("Option ", option ,", ", variant ,". variant",sep=""), ylab = "Variance", xlab = "lambda")
		#lines(lambda, 2*b/(b-1), col=4, lwd = 1)
		plot(lambda,  2*b/(b-1), type="l", lwd = 1, main = paste("Option ", option ,", ", variant ,". variant",sep=""), ylab = "Degree of freedom", xlab = "lambda")
		
		return(dof)
	}
	data_pos = c(4, 9, 14, 19, 24)
	dof = diag(sapply(data_pos, plot_ab))

	png(filename=paste("plot_tstat",option,".png",sep=""), width=1000, height=2500, pointsize=12)
	par(mfrow=c(5,2))
	# compare the t-statistic

	# defree of freedom, assumed
	#df 	= L + L - 2
	yy	= (1:(num*N))/(num*N)
	
	plot_tstat <- function(data_pos) {
		variant = (data_pos %/% 5) + 1	
		# turn list of t-statistic to matrix
		xts	= do.call("cbind", data[data_pos,])
		xx  = sort(as.vector(xts)) # for all genes, or xts[1,] for just the 1.gene	
		plot(xx, yy, type="l", lwd = 1, main = paste("Option ", option ,", ", variant ,". variant",sep=""), xlab = "");
		# variance of all t
		var_t = var(xx)
		df = 2*var_t/(var_t - 1)
		lines(xx, pt(xx, df), col=4, lwd=1)
		lines(xx, pt(xx, dof[variant]), col=2, lwd=1)
		plot(xx, yy-pt(xx, df), type="l", col=4, lwd=1, main = paste("Option ", option ,", ", variant ,". variant",sep=""), xlab = "");
		lines(xx, yy-pt(xx, dof[variant]), col=2, lwd=1)		
	}
	data_pos = c(3, 8, 13, 18, 23)
	sapply(data_pos, plot_tstat)
		
	png(filename=paste("hist_lambda",option,".png",sep=""), width=2000, height=1000, pointsize=12)
	# histogram for compare the lambda
	find_minmax_lambda <- function(data_pos) {
		c_max = c()
		c_min = c()
		for (i in 1:(length(data_pos))) {
			max = max(do.call("c", data[data_pos[i],]))
			min = min(do.call("c", data[data_pos[i],]))
			c_max = c(c_max, max)
			c_min = c(c_min, min)
		}
		x_max = max(c_max)
		x_min = min(c_min)
		return (list(x_max=x_max, x_min=x_min))
	}
	plot_lambda <- function(data_pos) {
		x_max = find_minmax_lambda(data_pos)$x_max
		x_min = find_minmax_lambda(data_pos)$x_min
		plot(NULL, NULL, xlim=c(x_min,x_max), ylim=c(0,1), main = paste("Option ",option,sep=""), xlab = "score", ylab = "Density")
		lg = c()
		for (i in 1:(length(data_pos))) {
			# calc which variant from the position of data in list
			variant = data_pos[i] - 25
			h 			= hist(do.call("c", data[data_pos[i],]), breaks = 100, plot = FALSE )
			h$density 	= h$counts/sum(h$counts)
			h$counts    = cumsum(h$counts)
			h$density   = cumsum(h$density)

			lines(h$mid, h$density, col=cols[i], lwd = 2)
			lg = c(lg, paste(variant ,". variant"))
		}

		legend(x="bottomright",inset=0.01,legend=lg,col=cols,lwd=3)
	}
	
	data_pos = c(26, 27, 28, 29, 30)
	plot_lambda(data_pos)
	
	png(filename=paste("splot_lambda_dof",option,".png",sep=""), width=1000, height=1500, pointsize=12)
	par(mfrow=c(3,2))
	# scatterplot lambda versus dof	
	splot_lambda_dof <- function(data_pos) {
		for (i in 1:(length(data_pos))) {
			# calc which variant from the position of data in list
			variant = data_pos[i] - 25
			# calc the position of the corresponded t-statistic
			pos_t = variant * 5 - 2
			lambda = do.call("c", data[data_pos[i],])
			# calc the degree of freedom from var
			dof = sapply(data[pos_t,], var)
			dof = 2 * dof / (dof - 1)
			if (length(dof) < length(lambda)) {
				dof = rep(dof, each = 2)
			}
			plot(lambda, dof, main = paste("Option ", option ,", ", variant ,". variant",sep=""), xlab = "lambda", ylab = "Dof")
		}
	}
#	data_pos = c(26, 27, 28, 29, 30)
	splot_lambda_dof(data_pos)

}

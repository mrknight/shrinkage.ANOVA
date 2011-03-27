
eval_data <- function(data, K, option) {

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
# load data
data_dir = "/Users/knight/dev/data/shrinkage_var/"
cols = rainbow(3*2)
simu_option1 = c("I", "II", "III-p001", "IV-p001")
simu_option2 = c("III-p002", "III-p005", "III-p01", "III-p02", "III-p05")
simu_option3 = c("IV-p002", "IV-p005", "IV-p01", "IV-p02", "IV-p05")

simu_option = c(simu_option1, simu_option2, simu_option3)


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
	cat(x_max, " ", x_min,"\n")
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
	cat(x_max, " ", x_min,"\n")

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


#data_pos = c(1, 6, 16)
#plot_score1(data_pos, info)

#simu_name = simu_option[10]

#load(paste(data_dir,"simu",simu_name,".rdata",sep=""))
#data = xd

#data_pos = c(2, 7, 12, 17, 22)
#plot_score_t(xd, data_pos, info = simu_name)	

plot_all_score <- function(simu_name) {
	load(paste(data_dir,"simu",simu_name,".rdata",sep=""))

	data_pos = c(1, 6, 16)
	plot_score(xd, data_pos, info = simu_name)	
}

plot_all_score_t <- function(simu_name) {
	load(paste(data_dir,"simu",simu_name,".rdata",sep=""))

	data_pos = c(2, 7, 12, 17, 22)
	plot_score_t(xd, data_pos, info = simu_name)	
}

#png("scoret_all.pdf")
#png(filename="scoret_all.png", width=1400, height=3500, pointsize=15)
#par(mfrow=c(7,2))

#sapply(simu_option, plot_all_score_t)

png(filename="score_all.png", width=1400, height=3500, pointsize=15)
par(mfrow=c(7,2))
sapply(simu_option, plot_all_score)

dof = rep(1, 5)	

num		= 500 # number of simulation
P		= 2000 # number of genes
N		= 3 # sample size
# defree of freedom, assumed
df2 	= N + N - 2
yy		= (1:(num*P))/(num*P)

eval_data <- function(data, option) {
	plot_ab <- function(data_pos) {
		# calc which variant from the position of data in list
		variant = (data_pos %/% 5) + 1

		# compute the mean of all a and b from all simulation
		a	= apply(do.call("cbind", data[data_pos,]), 1, mean)
		b	= apply(do.call("cbind", data[(data_pos+1),]), 1, mean)
	  	a[a < .Machine$double.eps] = 0
		# compute the defree of freedom from b
		dof[variant] = (mean(b)*2) / (mean(b) - 1)
		return(dof)
	}

	data_pos = c(4, 9, 14, 19, 24)
	dof = diag(sapply(data_pos, plot_ab))

	#pdf("a.pdf")
	png(filename=paste("plot_tstat",option,".png",sep=""), width=1700, height=2400, pointsize=15)
	par(mfrow=c(5,2))
	
	plot_tstat <- function(data_pos) {
		variant = (data_pos %/% 5) + 1	
		# turn list of t-statistic to matrix
		xts	= do.call("cbind", data[data_pos,])
		xx  = sort(as.vector(xts)) # for all genes, or xts[1,] for just the 1.gene	
		plot(xx, yy, type="l", lwd = 1, main = paste("Option ", option ,", ", variant ,". case",sep=""), xlab = "t-statistic", ylab = "Density", col = 1, cex.main = 2, cex.axis = 2, cex.lab = 2);
		# variance of all t
		var_t = var(xx)
		df = 2*var_t/(var_t - 1)
		lines(xx, pt(xx, df), col=2, lwd=1)
		lines(xx, pt(xx, dof[variant]), col=3, lwd=1)
		lines(xx, pt(xx, df2), col=4, lwd=1)

		lg = c(	"Cdf of shrinkage-t",
				"Cdf of t with estimated dof from variance", 
				"Cdf of t with estimated dof from lambda", 
				"Cdf of t with theoretical dof")
		cl = c(1:4)
#		legend(x="bottomright",inset=0.001,legend=lg,col=cl,lwd=2, cex = 2)
		plot(xx, yy-pt(xx, df), type="l", col=2, lwd=1, main = paste("Magnification for Option ", option ,", ", variant ,". case",sep=""), xlab = "t-statistic", ylab = "Density", ylim = c(-.05, .05)
			, cex.main = 2, cex.axis = 2, cex.lab = 2)
		lines(xx, yy-pt(xx, dof[variant]), col=3, lwd=1)
		lines(xx, yy-pt(xx, df2), col=4, lwd=1)
		abline(h = 0, col=1, lwd=2)

		legend(x="bottomright",inset=0.001,legend=lg,col=cl,lwd=5, cex = 2)

	#	ks = ks.test(yy, pt(xx, df))$statistic
	#	ks2 = ks.test(yy, pt(xx, df2))$statistic
	#	ks_dof = ks.test(yy, pt(xx, dof[variant]))$statistic
	#	cat(ks," ",ks2," ",ks_dof,"\n")
	}

	data_pos = c(3, 8, 13, 18, 23)
	sapply(data_pos, plot_tstat)

}

load_and_eval <- function(simu_name) {
	# loading data
	load(paste(data_dir,"simu",simu_name,".rdata",sep=""))

	eval_data(xd, simu_name)
}

#load_and_eval(simu_name)
sapply(simu_option1, load_and_eval)
#sapply(simu_option2, load_and_eval)
#sapply(simu_option3, load_and_eval)

dev.off()

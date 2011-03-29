source("/Users/knight/dev/src/shrinkage.ANOVA/lib.generate.R")

K = 2
p = 0
load(paste(dir.data,	"/K",K,"_p",p,".rdata", sep=""))

data = ed

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

dof 		= rep(0, 10)

data_pos 	= c(4, 9, 14, 19, 24)
dof 		= sapply(data_pos, plot_ab)


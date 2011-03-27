#################################
#
# \brief	the simulation for examine the performance of shrinkage ANOVA statistic
#
#################################

source("/Users/knight/dev/src/shrinkage.ANOVA/lib.generate.R")

m = matrix(1:15, 5, 3)
mm = c(1:5)
# testing 2
num 	= 10
K		= K_[3]

test2 <- function(l) {
	# generate samples
	X 	= sapply(1:P, generate_data_II, K, p = 0)
	
	f = eval.statistic(X, K, option = "")
	return (f)
}

x 		= sapply(1:num, test2)
xx 		= sort(as.vector(x))

z 		= (1:(num*P))/(num*P)
plot(xx, z, type = "l")
lines(xx, pf( xx, df1 = (K - 1), df2 = (N*K - K) ), col=3, lwd = 2)


# start simulation
#source("/Users/knight/dev/src/shrinkage_var/simu_shrinkage_t.R")
#data_dir = "/Users/knight/dev/data/shrinkage_var/"

#load_and_eval <- function(simu_name) {
	# loading data
#	load(paste(data_dir,"simu",simu_name,".rdata",sep=""))

#	eval_data(xd, simu_name)
#}

#sapply(simu_option1, load_and_eval)
#sapply(simu_option2, load_and_eval)
#sapply(simu_option3, load_and_eval)

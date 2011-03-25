#################################
#
# \brief	the simulation for examine the performance of shrinkage ANOVA statistic
#
#################################

source("/Users/knight/dev/src/shrinkage.ANOVA/lib.generate.R")

# testing 2
num 	= 20
K		= K_[3]

test2 <- function(l) {
	# generate samples
	X 	= sapply(1:P, generate_data_II, K, p = 0.001)
	
	f = eval.statistic(X, K, option = "")
	return (f)
}

x 		= sapply(1:num, test2)
xx 		= sort(as.vector(x))

z 		= (1:(num*P))/(num*P)
plot(xx, z, type = "l")
lines(xx, pf( xx, df1 = (K - 1), df2 = (N*K - K) ), col=3, lwd = 2)

simu_option1 = c("I", "II", "III-p001", "IV-p001")
simu_option2 = c("III-p002", "III-p005", "III-p01", "III-p02", "III-p05")
simu_option3 = c("IV-p002", "IV-p005", "IV-p01", "IV-p02", "IV-p05")

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

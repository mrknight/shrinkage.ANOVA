dir.source	= "/Users/knight/dev/src/shrinkage.ANOVA"
source(paste(dir.source, "/lib.generate.R"		, sep =""))

cols = rainbow(3*2)

#K	= K_[1]
#p	= p_[1]
#plot_all_score_t(K, p)
#data_dir = "/Users/knight/dev/data/shrinkage_var/"
#simu_option1 = c("I", "II", "III-p001", "IV-p001")
#simu_name = simu_option1[2]

#load(paste(data_dir,"simu",simu_name,".rdata",sep=""))
#data = xd

for (i in 1:length(K_)) {
	K	= K_[i]
	
#	png(filename=paste(dir.output,"/K",K,"scoret_all.png", sep=""), width=1400, height=2000, pointsize=15)
#	par(mfrow=c(4,2))	
	for (j in 1:length(p_)) {
		p	= p_[j]
#		plot_all_score_t(K, p)		
	}
	
#	png(filename=paste(dir.output,"/K",K,"score_all.png", sep=""), width=1400, height=2000, pointsize=15)
#	par(mfrow=c(4,2))	
	for (j in 1:length(p_)) {
		p	= p_[j]
#		plot_all_score(K, p)	
	}
	
	for (j in 1:length(p_)) {
		# loading data
		p 	= p_[j]
		load(paste(dir.data,	"/K",K,"_p",p,".rdata", sep=""))
		output.F.data(ed, K, p)
	}
}
		
#png(filename="score_all.png", width=1400, height=3500, pointsize=15)
#par(mfrow=c(7,2))
#sapply(simu_option, plot_all_score)

#dof = rep(1, 5)	

#num		= 500 # number of simulation
#P		= 2000 # number of genes
#N		= 3 # sample size
# defree of freedom, assumed
#df2 	= N + N - 2
#yy		= (1:(num*P))/(num*P)
#load_and_eval(simu_name)
#sapply(simu_option1, load_and_eval)

dev.off()
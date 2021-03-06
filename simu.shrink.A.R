dir.source	= "/Users/knight/dev/src/shrinkage.ANOVA"
source(paste(dir.source, "/lib.generate.R"		, sep =""))

#eval.F.statistic(1, 5, 0)
for (i in 1:length(K_)) {
	for (j in 1:length(p_)) {
		K	= K_[i]
		p	= p_[j]
		start 		= proc.time() # time start
				
		cat("Generate data for K =",K,"with p =",p,"outliers...\n")
		# generate the data and evaluate on the fly
		ed	= sapply(1:num, eval.F.statistic, K, p)

		# save evaluated data to file
		save(ed, file = paste(dir.data,	"/K",K,"_p",p,".rdata", sep=""))
		
		time 		= proc.time() - start # time stop
		cat("Done in",time[3],"s...\n")
	}
}

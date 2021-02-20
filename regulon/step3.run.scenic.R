	.libPaths("/PRJ/software/R_library351")

	load("int/expression_MF.RData")
cat("###############################################################################\n")
	#print(sessionInfo())
cat("###############################################################################\n")
cat(date(),": about to run SCENIC\n")


	library(SCENIC)
	scenicOptions <- readRDS("int/scenicOptions.Rds")
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$nCores <- 6
	scenicOptions@settings$seed <- 123
	#scenicOptions@dbs@500bp<- "hg19-500bp-upstream-7species.mc9nr.feather"
	#scenicOptions@dbs@10kb<- "hg19-tss-centered-10kb-7species.mc9nr.feather"


cat("###############################################################################\n")
cat(date(), ": about to run runSCENIC_1_coexNetwork2modules\n")
	
	runSCENIC_1_coexNetwork2modules(scenicOptions)



cat("###############################################################################\n")
cat(date(),": about to run runSCENIC_2_createRegulons\n")


	runSCENIC_2_createRegulons(scenicOptions)



cat("###############################################################################\n")
cat(date(),": about to run runSCENIC_3_scoreCells\n")
	
	runSCENIC_3_scoreCells(scenicOptions, exprMat_L)

	cat("finished\n")

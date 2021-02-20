	load("int/expression_MF.RData")
	.libPaths("/home/groups/kathryncheah/software/R_library351")
cat("###############################################################################\n")
cat(date(), ": about to run corrMat\n")

	library(SCENIC)
	org="mgi" # or hgnc, or dmel
	#dbDir="../scenic_databases" # RcisTarget databases location
	#myDatasetTitle="SCENIC example on Mouse brain" # choose a name for your analysis
	myDatasetTitle<-"SCENIC on irene" # choose a name for your analysis
	dbDir<-"/home/groups/kathryncheah/pkchen/hisat2/outputs-grch38-single-cell/analysis-2531/scenic_databases/" # RcisTarget databases location
	scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=10) 

	#scenicOptions@settings$nCores <- 6
	scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
	scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

	saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
	#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
	 #                          minCountsPerGene=3*.01*ncol(exprMat),
	  #                         minSamples=ncol(exprMat)*.01)

	exprMat2<-t(exprMat_filtered)
	exprMat2[exprMat2==0]<-NA

	corrMat <- cor(exprMat2, method="spearman",use="pairwise.complete.obs")
	saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

cat("###############################################################################\n")
	cat(date(), ": about to run runGenie3\n")
	runGenie3(exprMat_L, scenicOptions)


cat("finished\n")

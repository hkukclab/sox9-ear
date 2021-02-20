.libPaths("/home/groups/kathryncheah/software/R_library351")

	if(1){
		dir.create("int")
		dir.create("output")
	}
	options(stringsAsFactors=FALSE)


	library(SCENIC)
	org<-"mgi" # or hgnc, or dmel
	dbDir<-"/home/groups/kathryncheah/pkchen/hisat2/outputs-grch38-single-cell/analysis-2531/scenic_databases/" # RcisTarget databases location
	myDatasetTitle<-"SCENIC on irene" # choose a name for your analysis
	scenicOptions <- initializeScenic(org=org, dbDir=dbDir,
		 datasetTitle=myDatasetTitle, nCores=12) 
cat("###############################################################################\n")
	load("../step1.cell.infor.246.slim.RData")
	load("../cell.pheno.in.mysql.246.aug13.2019.RData")
	#tmp<-read.table("../siglist/all.DEGs")
	#SigAll<-as.character(tmp[,1])

	#ind.sc<-which(cell.pheno$is.bulk=="FALSE")
#	flag1<-cell.pheno$is.bulk=="FALSE"
#	flag2<-lab1$newLab1!="others"

	indMF<-colSums(dat.counts.mRNA>0)>1e3
	cat("length MF:",table(indMF),"\n")
	#which(flag1&flag2)
	geneUNQ<-unique(annot.irene[,4])
	#geneUNQ<-intersect(SigAll,unique(annot.2531[,4]))
#	geneUNQ<-intersect(SigAll,unique(annot.irene[,4]))
	indUNQ<-match(geneUNQ,annot.irene[,4])
	exprMat<-dat.counts.mRNA[indUNQ,indMF]
	rownames(exprMat)<-annot.irene[indUNQ,4]

#	print(table(lab1$patient[indMF],lab1$newLab1[indMF]))
	print(head(exprMat[,1:4]))
	
cat("###############################################################################\n")
	cellInfo <- droplevels(dat.patient[indMF,])
	head(cellInfo)
	rownames(cellInfo)<-colnames(exprMat)
	#table(cellInfo$)


	
	#genesKept<- rowSums(exprMat>0)>100
	genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=1,
                           minSamples=15)

	exprMat_filtered <- exprMat[genesKept,]
	cat("dim of exprMat_filtered:",dim(exprMat_filtered),"\n")

	exprMat_L <- log2(exprMat_filtered+1) 
#	print(as.matrix(table(cell.pheno$n.patient[indMF])))



	save(exprMat,cellInfo,
		exprMat_filtered,
		exprMat_L,
		file="int/expression_MF.RData")


cat("###############################################################################\n")
	saveRDS(cellInfo, file="int/cellInfo.Rds")
	colVars <- list(CellType=setNames(c("forestgreen", "darkorange", "magenta4", "hotpink", "red3", "skyblue", "darkblue","darkgrey"), 
	    c("MRC_CD", "MRC_WT", "Prolif_CD", "Prolif_WT", "RRC_CD", "RRC_WT", "","---")))
	saveRDS(colVars, file="int/colVars.Rds")

	cat("Step-1 finished \n")

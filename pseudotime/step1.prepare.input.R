
load("../resources/step1.cell.infor.246.slim.RData")
load("../resources/cell.pheno.in.mysql.246.oct9.2019.RData")



library(monocle)
cellnames<-colnames(dat.counts.mRNA)

table(cellnames==dat.patient$sample.ids)
#write.table(data.frame(genes[,1:8],fpkmMat.new),file="Gene.annot.fpkm.469.txt",quote=F,row.names=F,sep="\t")
#write.table(data.frame(annot.gtf[,c(1,2,4)],dat.counts.new),file="Gene.annot.469.txt",quote=F,row.names=F,sep="\t")
#write.table(cell.annot.new,file="cell.annot.469.txt",quote=F,row.names=F,sep="\t")

#############################

FENO<-droplevels(data.frame(cellname=dat.patient$sample.ids,
	dat.patient[,-2],
	rate=as.numeric(gsub("%","",alignment[,12]))))
RAWCOUNTS<-dat.counts.mRNA
colnames(RAWCOUNTS)<-rownames(FENO)

GENES<-data.frame(annot.irene,gene_short_name=annot.irene[,4])
rownames(GENES)<-annot.irene[,1]

pd <- new("AnnotatedDataFrame", data = FENO)
fd <- new("AnnotatedDataFrame", data = GENES)

HSMM <- newCellDataSet(as.matrix(RAWCOUNTS),
		 phenoData = pd, featureData = fd,
	    expressionFamily=negbinomial.size())

save(HSMM,alignment,file="HSMM.htseq.246.RData")


#############################
ind.bulk<-which(dat.patient$cellsource=="BULK")

FENO<-droplevels(data.frame(cellname=cellnames,
	dat.patient,
	rate=as.numeric(gsub("%","",alignment[,12]))))[-ind.bulk,]
RAWCOUNTS<-dat.counts.mRNA[,-ind.bulk]
colnames(RAWCOUNTS)<-rownames(FENO)

GENES<-data.frame(annot.kwok,gene_short_name=annot.kwok[,4])
rownames(GENES)<-annot.kwok[,1]

pd <- new("AnnotatedDataFrame", data = FENO)
fd <- new("AnnotatedDataFrame", data = GENES)

HSMM <- newCellDataSet(as.matrix(RAWCOUNTS),
		 phenoData = pd, featureData = fd,
	    expressionFamily=negbinomial.size())

save(HSMM,alignment,file="HSMM.htseq.443.RData")



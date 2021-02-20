
load("../resources/step1.cell.infor.246.slim.RData")
load("../resources/cell.pheno.in.mysql.246.aug13.2019.RData")

table(colnames(fpkmMat)==dat.patient$sample.ids)

ind.cells<-which(colSums(fpkmMat>0)>1e3)
GENES<-genes[,1:7]

#################################
#################################
library(MAST)
source("../resources/DEGs.single.cell.MAST.R")


L.res.all<-list()
L.IND<-L.res.all


	cell_label<-dat.patient$tissue[ind.cells]
	print(as.matrix(table(cell_label)))


	scMAT<-fpkmMat[,ind.cells]
	FENO<-data.frame(colnames(scMAT),cell_label)

	cnt<-cnt+1
	MYDEG("results-allCD-vs-allWT",
		"WT",
		"CD",
		"allWT",
		"allCD",
		1,1,SPECIES="mouse")



save(L.res.all,L.IND,
	levpatients,
	file="step1.allCD-vs-allWT.aug2019.RData")
#################################
#################################

library(gplots)

source("../resources/addTFCD.R")

########################################
fn1<-"step1.allCD-vs-allWT.aug2019.RData"

library(grid)
library(lattice)
library(gridExtra)
library(stringr)

tsne246<-read.table("../resources/tsne-246/pca.2pc.thismarkers.for.irene_CD_WT_3.txt")
load("../resources/cell.pheno.in.mysql.246.aug13.2019.RData")

########################################
pdf("volcano.allCD-vs-allWT.feb2021.pdf",height=6,width=14)
	load(fn1)
	for(i in 1:length(L.res.all)){
		DEGi<-L.res.all[[i]]

		layout(t(matrix(seq(3))),widths=c(3,5.5,5.5))
		thisType<-as.character(DEGi$group2[1])
		otherType<-as.character(DEGi$group1[1])

		par(mar=c(12,4,12,0))
		indThis<-which(dat.patient$tissue==thisType)
		indOthers<-which(dat.patient$tissue==otherType)
		indi<-c(indThis,indOthers)

		tsnei<-tsne246[match( as.character(dat.patient$sample.ids),rownames(tsne246)),]

		MEMBi<-rep(c(thisType,otherType),c(length(indThis),length(indOthers)))
		
		
		TABi<-table(MEMBi)
		plot(tsnei,pch=c(18,16)[as.factor(MEMBi)],col="grey",cex=1,
			main=paste0("between ",TABi[1]," ",names(TABi)[1],
			" and ",TABi[2]," ",names(TABi)[2]))
		points(tsnei[c(indThis,indOthers),],cex=3,
			pch=c(18,16)[as.factor(MEMBi)],
			col=c("#F8766D","#619CFF")[as.factor(MEMBi)])
		legend("bottomleft",legend=paste0(names(TABi)," ",c("(RED)","(BLUE)"),": ",TABi),bty='n')



		par(mar=c(5,4,4,2))

		irrGene<-!grepl("^RP|\\.|-",DEGi$gene_short_name)

		rowMax<-apply(DEGi[,7:8],1,max)
		gap<-DEGi$gap
		pval<-DEGi$fdr
		genei<-Add_TFCD(as.character(DEGi$gene_short_name))

		indUp<-which(gap>0.2  & DEGi$fdr<0.05&rowMax>1&irrGene)
		indDown<-which(gap< -0.2 & DEGi$fdr<0.05&rowMax>1&irrGene)
		indSig<-c(indUp,indDown)
		if(length(indSig)==0){
			YLIM<-c(0,0)
		}else{
			YLIM<-c(0,max(-log10(pval[c(indUp,indDown)])))
		}
		plot(gap,-log10(pval),pch=16,cex=1/2,col="grey",
			xlab="",
			ylab="",
			ylim=YLIM*1.2,type='n',
			main="all WT (blue) vs all CD (red)")
		points(gap[-indSig],-log10(pval)[-indSig],pch=16,cex=1,col="grey",
			xlab="",
			ylab="",
			ylim=YLIM*1.2,type='n',
			main="")


		abline(v=0,lty=2)
		abline(v=c(-0.2,0.2),lty=3)
		abline(h=-log10(0.05),lty=3)
		mtext(paste0("RED",
			" (", DEGi$group2[1], DEGi$n2[1]," cells)  vs.  ",
			"BLUE",
			" (", DEGi$group1[1], DEGi$n1[1]," cells)"),
			cex=3/4)

		if(length(indUp)>0){
			points(gap[indUp],-log10(pval)[indUp],pch=16,cex=abs(DEGi$LOGFC)[indUp],col="#F8766D")
			text(gap[indUp],-log10(pval)[indUp], 
				genei[indUp],pos=rep(c(2,4),100),cex=1,offset=0.2,xpd=T)

		
		}
		if(length(indDown)>0){
			indTF<-grep("\\*|#",genei[indDown])
			points(gap[indDown],-log10(pval)[indDown],pch=16,cex=abs(DEGi$LOGFC)[indDown],col="#619CFF")
			text(gap[indDown],-log10(pval)[indDown], 
				genei[indDown],pos=rep(c(2,4),100),cex=1,offset=0.2,xpd=T)
		}

		legend("top",pch=16,col=c(2,4),
			legend=c(paste0("Gap>0.2: specifically high to red cells (n=",length(indUp),")"),
				paste0("Gap< -0.2: specifically low to red (up in blue) cells (n=",length(indDown),")")),
			cex=2/3,bty='n')
		legend("top",pch=16,pt.cex=seq(4)/2,inset=c(0,0.1,0,0),
			title="||log2(FC)||",
			legend=seq(4)/2,bty='n',ncol=4)

		##################
		genei[genei=="FBN1*"]<-"FBN1"

		gi<-unique(genei[indDown]);
		TFdown<-gsub("\\*","",paste0(gi[grep("\\*",gi)],collapse=", "))

		CDdown<-gsub("#","",paste0(gi[grep("#",gi)],collapse=", "))

		othergene<-setdiff(gi[!grepl("#",gi)&!grepl("\\*",gi)],c("XIST","TSIX"))
		otherDown<-paste0(othergene,collapse=", ")

		##################

		gi<-unique(genei[indUp]);
		TFup<-gsub("\\*","",paste0(gi[grep("\\*",gi)],collapse=", "))

		CDup<-gsub("#","",paste0(gi[grep("#",gi)],collapse=", "))

		othergene<-setdiff(gi[!grepl("#",gi)&!grepl("\\*",gi)],c("XIST","TSIX"))
		otherUp<-paste0(othergene,collapse=", ")


		##################
		mytheme <- gridExtra::ttheme_default(base_size =8,padding = unit(c(1, 1), "mm"))

		WIDTH<-40
		mat1<-rbind(c(str_wrap(TFdown,WIDTH),str_wrap(TFup,WIDTH)),
			c(str_wrap(CDdown,WIDTH),str_wrap(CDup,WIDTH)),
			c(str_wrap(otherDown,WIDTH),str_wrap(otherUp,WIDTH)))
		rownames(mat1)<-c("TF","CD","others")
		colnames(mat1)<-c("up in blue","up in red")
		pushViewport(viewport(x=0.8,y=0.5))
		grid.table(mat1,theme=mytheme )
		popViewport()
	}
dev.off()

#################################








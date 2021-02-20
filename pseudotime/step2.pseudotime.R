load("HSMM.htseq.246.RData")

library(monocle2)

HSMM2 <- detectGenes(HSMM, min_expr = 1)
valid_cells <- rownames(subset(pData(HSMM2 ), 
	num_genes_expressed>1e3&rate>50))

HSMM3 <- HSMM2[,valid_cells]

#############################

expressed_genes <- row.names(subset(fData(HSMM3 ), num_cells_expressed >= 10))

HSMM4 <- estimateSizeFactors(HSMM3)
HSMM5 <- estimateDispersions(HSMM4)

#############################
 set.seed(0)
#############################
disp_table <- dispersionTable(HSMM5)

ordering_genes <- droplevels(subset(disp_table,
	mean_expression >= 1 &
	dispersion_empirical >= 2 * dispersion_fit)$gene_id)


HSMM6 <- setOrderingFilter(HSMM5, ordering_genes)
#############################



#############################################
HSMM7<- reduceDimension(HSMM6, max_components=3)
HSMM8<- orderCells(HSMM7, reverse=FALSE)
#############################################
resplot<-plot_cell_trajectory(HSMM8,color_by="ireneLabel")

plot_cell_trajectory(HSMM8,color_by="ireneLabel",cell_size = 8)

pdf("ear.246.pseudotime.pdf")
	plot_cell_trajectory(HSMM8,color_by="Pseudotime",cell_size = 8)
	plot_cell_trajectory(HSMM8,color_by="ireneLabel",cell_size = 8)
	plot_cell_trajectory(HSMM8,color_by="IreneCond2",cell_size = 8)
	plot_cell_trajectory(HSMM8,color_by="cellsource",cell_size = 8)
	plot_cell_trajectory(HSMM8,color_by="State",cell_size = 8)
	plot_cell_trajectory(HSMM8,color_by="num_genes_expressed",cell_size = 8)

dev.off()








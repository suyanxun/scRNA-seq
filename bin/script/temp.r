##单独运行
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

#features<-c("CEACAM1","CEACAM6","CEACAM3","MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","THY1","ITGA6","CD38","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CCR7","CD4","FOXP3","SELL","S100A4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","LILRA4","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")

HSC = c("NPTX2", "H1F0", "EMP1", "MEIS1", "CALCRL", "TPSD1", "TPT1", "CRHBP", "CLNK", "TSC22D1", "DST", "NRIP1", "ABCB1", "GABRA4", "ZBTB20", "ABCA9", "TPSB2", "KMT2A", "FAM30A", "MEF2C", "TMEM74", "PDZRN4", "ST3GAL1", "XIRP2", "RBPMS", "TMEM25", "C20orf203", "GNG11", "SLC6A13", "HOPX")
Progenitor = c("CDK6", "HSP90AB1", "SPINK2", "EEF1B2", "PCNP", "TAPT1-AS1", "HINT1", "LRRC75A-AS1", "DSE", "PEBP1", "LOC107984974", "H2AFY", "EEF1A1", "SMIM24", "PSME1", "SOX4", "LINC01623", "EEF1G", "EBPL", "EIF4B", "PARP1", "MEST", "TMEM70", "TFDP2", "ATP5G2", "NAP1L1", "MSI2", "TPM4", "SPN", "SELL")
GMP = c("PRTN3", "MPO", "CALR", "CLEC5A", "ELANE", "POU4F1", "TRH", "TSPOAP1", "CEBPE", "LINC01835", "NUCB2", "CSF3R", "RUNX1T1", "CD38", "PLPPR3", "IGFBP2", "PRRT4", "SNHG5", "FABP5", "LOC100419170", "CLEC11A", "SERPINB1", "AZU1", "FBN2", "HNRNPDL", "HSPB1", "RNA5-8S", "THSD7A", "C12orf57", "FGFR1")
Promono = c("DEFB1", "RNASE2", "MS4A3", "SERPINB10", "SESN3", "ZFR", "MRPL33", "CTSG", "SLC44A1", "SLPI", "FUT4", "SRGN", "CD70", "PRLR", "PLD3", "LPL", "RETN", "TP53INP2", "HSPA5", "RNASE3", "CCL23", "EMB", "ATP8B4", "CLU", "FAM107B", "KBTBD11", "CSTA", "ANKRD28", "PIWIL4", "RNVU1-6")
Monocyte = c("FCN1", "S100A12", "MAFB", "VCAN", "S100A9", "PLBD1", "SERPINA1", "BCL2A1", "THBS1", "PSAP", "S100A8", "FPR1", "C5AR1", "CD14", "NAMPT", "VNN2", "CTSS", "DUSP1", "CEBPB", "CR1", "NFKBIA", "SLC11A1", "LILRB3", "BCL6", "CYP1B1", "TNFAIP2", "MS4A10", "AQP9", "TLR4", "APOBEC3A")
cDC =  c("MRC1", "HLA-DRB5", "CST3", "SAMHD1", "NAPSB", "FCER1A", "HLA-DRB1", "JAML", "PKIB", "HLA-DRA", "HLA-DRB6", "CPVL", "HLA-DPB1", "HLA-DQA1", "HLA-DPA1", "CLEC4A", "TMSB10", "CAP1", "HLA-DQB1", "CRIP1", "CLEC10A", "GPX1", "ITGB7", "HLA-DQB2", "DBI", "FTH1P3", "ACTB", "HLA-DQA2", "S100B", "ALDH2")

cell_normal = c( "CSF1", "MNDA", "FCER1A", "FCER1G", "CD34", "LYST", "MEIS1", "CD8A", "MME", "CD38", "PROM1", "KIT", "JCHAIN", "BANK1", "GYPA", "IL7R", "GZMK", "TCF7", "EGR1", "MZB1", "CD14", "EBF1", "CD24", "CEBPD", "PAX5", "FCN1", "EGFL7", "HBB", "HBD", "MS4A1", "NCAM1", "CD3D", "CD3G", "CLEC4C", "CLEC4A", "KLRB1", "KLRD1", "LYZ", "CTSG", "GZMB", "IL32", "CD19", "IRF8", "CLEC10A", "CCL5", "MSI2", "MPO", "TCF4", "AZU1", "ELANE", "PTPRS", "CD79A", "C5AR1", "VPREB1", "IGLL5" )

cell_normal = c( "MEIS1", "EGR1", "MSI2", "CD38", "CD34", "PROM1", "EGFL7", "MPO", "ELANE", "CTSG", "AZU1", "LYST", "LYZ", "CEBPD", "MNDA", "FCER1G", "FCN1", "CD14", "C5AR1", "CLEC4A", "CLEC10A", "FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4", "CSF1", "KIT", "HBB", "HBD", "GYPA", "CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19", "CD79A", "MS4A1", "BANK1", "MZB1", "IGLL5", "JCHAIN", "CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5", "GZMK", "CD8A", "KLRB1", "KLRD1", "GZMB", "NCAM1")

flowCytometryFeatures<-c("PTPRC","CD19","MME","MS4A1","CD27","CR2","CD38","IL3RA","ITGAX","HLA-DRA","HLA-DRB1","CD3D","CD3E","CD3G","CD14","FCGR3A","THBD","CD1C","CD34","KIT","CD33","ANPEP","ITGAM","FCGR1A","NCAM1","CD7","CD8A","CD8B","CCR7","CD69","CD4","ITGAE")

dims=c(10,15,16,17,18,19,20,21)
#dims=c(20)
res=c(0.5,0.8,1.0,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8)


width = 6
height = 6
sample_temp = "all"
sample = sample_temp

features = cell_normal
type = "cell_normal"

indir = paste( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/", sample, sep = "")
outdir = paste( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/", sample, "/", type, sep = "")

#indir = paste( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/JCI/Result/", sample, sep = "")
#outdir = paste( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/JCI/Result/", sample, sep = "")

sample_name = 1:15

for ( dim in dims ){
	sample = paste( sample_temp, "_Dim", dim, sep = "" )
	seurat_object_file = paste( indir, "/", sample, "_cluster.rds", sep = "")
	print( seurat_object_file )
	seurat_object = readRDS( seurat_object_file )
	p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "stim", pt.size = 0.1 )
	p2 <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.1 )
	plots = plot_grid(p1, p2)
	ggsave( paste( outdir, "/", sample, "_cluster_DimPlot.pdf" ,sep = "" ), plots, width = width, limitsize = FALSE , height = height)
	plots <- DimPlot(seurat_object, reduction = "umap", group.by = "stim", pt.size = 0.1, width = width, height = height  ) + NoLegend()
	ggsave( paste( outdir, "/", sample, "_cluster_DimPlot_NoLegend.pdf" ,sep = "" ), plots, width = width,limitsize = FALSE, height = height )
	plots = DimPlot(seurat_object, reduction = "umap", split.by = "stim" )
	ggsave( paste( outdir, "/", sample, "_sample_DimPlot.pdf" ,sep = "" ), plots, width = width*2,limitsize = FALSE, height = height*2 )
	nrow <- ceiling(length(res)/4)
	p<-1:length(res)
	p <- paste(rep("p",length(res)),p,sep = "")
	for ( i in 1:length(res) ){
		if (length(sample_name) == 1){
			col_name = paste( "RNA_snn_res.", res[i], sep = "" )
		}else{
			col_name = paste( "integrated_snn_res.", res[i], sep = "" )
		}
		
		temp <- DimPlot(object = seurat_object, reduction="umap", pt.size = 0.1, label = TRUE, label.size=6,group.by = col_name,order=c(levels(factor(seurat_object@meta.data[col_name],levels=c((dim(unique(seurat_object@meta.data[col_name]))[1]-1):0))))) +NoLegend() + labs(title=col_name)
		assign( p[i], temp )
	}
	text = paste( paste( "get(p[",1:length(res),"]) ", sep = "" ), collapse = ",")
	eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=4 )",sep = "" )))
	ggsave( paste( outdir, "/", sample, "_Different_Res_DimPlot.png",sep = "" ), plots, width = 4*width,height = height*nrow, limitsize = FALSE)
	features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
	plots = VlnPlot(seurat_object, features = features, ncol=4, pt.size = 0.1)
	ggsave( paste( outdir, "/", sample, "_biomarkers_VlnPlot.png" ,sep = "" ), plots, height = height*floor(length(features)/4),limitsize = FALSE, width = 4*width   )
	plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9" )
	ggsave( paste( outdir, "/", sample, "_biomarkers_FeaturePlot.png" ,sep = "" ), plots, height = height*floor(length(features)/4),limitsize = FALSE, width = 4*width  )
}

q()

library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

seurat_object = readRDS("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/P105SR_25SampleAlignmentSeurat.rds")

express_matrix = t(as.matrix(seurat_object[["RNA"]]@data))
gene_list = colnames(express_matrix)
cell_ID = rownames(express_matrix)
cell_type = as.character(seurat_object@active.ident[cell_ID])
cluster = as.character(seurat_object@meta.data$integrated_snn_res.3.8)
phase = as.character(seurat_object@meta.data$Phase)
sample = as.character( seurat_object@meta.data$stim )

write.table( express_matrix, "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/express_matrix.txt", sep="\t",quote = FALSE )
write.table( paste(gene_list, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/gene_list.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )
write.table( paste(cell_ID, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/cell_ID.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )
write.table( paste(cell_type, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/cell_type.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )
write.table( paste(cluster, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/cluster.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )
write.table( paste(phase, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/phase.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )
write.table( paste(sample, collapse="\t"), "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P105/SPRING/sample.txt", sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE  )


plots = DimPlot(seurat_object, reduction = "umap", split.by = "response",pt.size = 0.1) + NoLegend()
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_response_DimPlot.pdf" , plots, width = 30, limitsize = FALSE, height = 10 )

plots = DimPlot(seurat_object, reduction = "umap", split.by = "sample",pt.size = 0.1) + NoLegend()
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_sampleName_DimPlot.pdf" , plots, width = 30, limitsize = FALSE, height = 30 )

plots = DimPlot(seurat_object, reduction = "umap", split.by = "FAB_class",pt.size = 0.1) + NoLegend()
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_FAB_DimPlot.pdf" , plots, width = 30, limitsize = FALSE, height = 20 )


###############差异表达基因分析############################
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
seurat_object = readRDS("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds")


Idents(seurat_object) = seurat_object[["FAB_class"]]
nk.markers = FindAllMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, grouping.var = "stim")
write.table(nk.markers ,file = paste( opt$o, "/", file_name, "_all_marker.txt" ,sep = "" ), quote = FALSE,row.names = TRUE, sep = "\t")


write.table(nk.markers ,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/DE_gene/cellType_Marker.txt", quote = FALSE,row.names = TRUE, sep = "\t")




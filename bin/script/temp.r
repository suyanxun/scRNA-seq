##单独运行
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

features<-c("CEACAM1","CEACAM6","CEACAM3","MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","THY1","ITGA6","CD38","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CCR7","CD4","FOXP3","SELL","S100A4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","LILRA4","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")

#cell_normal = c( "CSF1", "MNDA", "FCER1A", "FCER1G", "CD34", "LYST", "MEIS1", "CD8A", "MME", "CD38", "PROM1", "KIT", "JCHAIN", "BANK1", "GYPA", "IL7R", "GZMK", "TCF7", "EGR1", "MZB1", "CD14", "EBF1", "CD24", "CEBPD", "PAX5", "FCN1", "EGFL7", "HBB", "HBD", "MS4A1", "NCAM1", "CD3D", "CD3G", "CLEC4C", "CLEC4A", "KLRB1", "KLRD1", "LYZ", "CTSG", "GZMB", "IL32", "CD19", "IRF8", "CLEC10A", "CCL5", "MSI2", "MPO", "TCF4", "AZU1", "ELANE", "PTPRS", "CD79A", "C5AR1", "VPREB1", "IGLL5" )

cell_normal = c( "MEIS1", "EGR1", "MSI2", "CD38", "CD34", "PROM1", "EGFL7", "MPO", "ELANE", "CTSG", "AZU1", "LYST", "LYZ", "CEBPD", "MNDA", "FCER1G", "FCN1", "CD14", "C5AR1", "CLEC4A", "CLEC10A", "FCER1A", "CLEC4C", "PTPRS", "IRF8", "TCF4", "CSF1", "KIT", "HBB", "HBD", "GYPA", "CD24", "EBF1", "MME", "VPREB1", "PAX5", "CD19", "CD79A", "MS4A1", "BANK1", "MZB1", "IGLL5", "JCHAIN", "CD3D", "CD3G", "IL32", "IL7R", "TCF7", "CCL5", "GZMK", "CD8A", "KLRB1", "KLRD1", "GZMB", "NCAM1")

# other_marker = c( "ITGA2B", "PF4", "VWF", "CA1", "HBB", "KLF1", "TFR2", "CCR2", "IRF8", "MPEG1", "ELANE", "MPO", "LYZ", "CSF1R", "CTSG", "PRTN3", "AZU1", "RGS1", "NPTX2", "DDIT4", "ID2", "DNTT", "RAG1", "RAG2", "CRHBP", "HLF", "DUSP1", "KLF1", "CA1", "ITGA2B", "PLEK", "CLC", "CPA3","HDC", "DNTT", "CD79A", "VPREB1", "IRF8", "SPIB", "IGKC" )


flowCytometryFeatures<-c("PTPRC","CD19","MME","MS4A1","CD27","CR2","CD38","IL3RA","ITGAX","HLA-DRA","HLA-DRB1","CD3D","CD3E","CD3G","CD14","FCGR3A","THBD","CD1C","CD34","KIT","CD33","ANPEP","ITGAM","FCGR1A","NCAM1","CD7","CD8A","CD8B","CCR7","CD69","CD4","ITGAE")

#MRD_Marker = c( "CD56", "CD13", "CD123", "CD4", "CD64", "CD133", "CD7", "CD41", "CD71", "NG2" )
MRD_Marker = c( "NCAM1", "ANPEP", "IL3RA", "CD4", "FCGR1A", "PROM1", "CD7", "ITGA2B", "TFRC", "CSPG4" )

#dims=c(10,15,16,17,18,19,20,21)
dims=c(20)
res=c(0.5,0.8,1.0,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8)


width = 6
height = 6
sample_temp = "all"
sample = sample_temp

features = MRD_Marker
type = "MRD_Marker"

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

##################
for (i in names(table(sample))){
	cells = names(sample[sample == i])
	plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",cells = c(sample == i) )
	ggsave( paste( outdir, "/", sample, "_biomarkers_FeaturePlot","_",i,.png" ,sep = "" ), plots, height = height*floor(length(features)/4),limitsize = FALSE, width = 4*width  )
}
##################
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


##############################################查看cell中的文章的数据情况###########################################
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
library(pheatmap)

HSC = c("HOPX", "TSC22D1", "FAM30A", "MEF2C", "TPSB2", "GNG11", "TPT1", "TPSD1", "DST", "NRIP1", "ABCB1", "ZBTB20", "KMT2A", "ST3GAL1", "RBPMS", "TMEM74", "H1F0", "EMP1", "MEIS1", "PDZRN4", "XIRP2", "TMEM25", "C20orf203", "SLC6A13", "NPTX2", "ABCA9", "GABRA4", "CALCRL", "CLNK", "CRHBP")

Progenitor = c("CDK6", "HSP90AB1", "EEF1B2", "PCNP", "HINT1", "LRRC75A-AS1", "PEBP1", "LOC107984974", "H2AFY", "EEF1A1", "SMIM24", "PSME1", "SOX4", "EEF1G", "EBPL", "EIF4B", "PARP1", "MEST", "TMEM70", "TFDP2", "ATP5G2", "NAP1L1", "MSI2", "TPM4", "SPN", "SPINK2", "SELL", "TAPT1-AS1", "DSE", "LINC01623")

GMP = c("MPO", "SERPINB1", "CLEC11A", "AZU1", "CALR", "NUCB2", "CSF3R", "ELANE", "TRH", "RUNX1T1", "CD38", "POU4F1", "CEBPE", "LINC01835", "SNHG5", "FABP5", "LOC100419170", "HNRNPDL", "HSPB1", "RNA5-8S", "C12orf57", "THSD7A", "PRTN3", "CLEC5A", "PLPPR3", "IGFBP2", "PRRT4", "FBN2", "TSPOAP1", "FGFR1")

Promono = c("SRGN", "CSTA", "MRPL33", "ANKRD28", "RNASE2", "PLD3", "RETN", "HSPA5", "EMB", "ZFR", "FAM107B", "MS4A3", "CTSG", "SLC44A1", "SLPI", "FUT4", "RNASE3", "CCL23", "ATP8B4", "CLU", "KBTBD11", "PIWIL4", "DEFB1", "SERPINB10", "SESN3", "CD70", "PRLR", "LPL", "TP53INP2", "RNVU1-6")

Monocyte = c("S100A9", "S100A8", "PSAP", "FCN1", "S100A12", "SERPINA1", "BCL2A1", "FPR1", "CTSS", "DUSP1", "CEBPB", "NFKBIA", "VCAN", "CD14", "NAMPT", "TNFAIP2", "SLC11A1", "LILRB3", "BCL6", "CYP1B1", "MS4A10", "AQP9", "TLR4", "MAFB", "PLBD1", "THBS1", "C5AR1", "VNN2", "CR1", "APOBEC3A")

cDC =  c("HLA-DQA1", "HLA-DRA", "HLA-DRB6", "HLA-DPB1", "HLA-DRB1", "DBI", "CRIP1", "HLA-DQB1", "HLA-DPA1", "CAP1", "TMSB10", "SAMHD1", "JAML", "CST3", "ACTB", "NAPSB",  "CLEC10A", "GPX1", "ITGB7",  "FTH1P3", "FCER1A",  "MRC1", "HLA-DRB5", "S100B", "CPVL", "ALDH2", "CLEC4A", "PKIB", "HLA-DQA2", "HLA-DQB2")
#########################
HSC = c("HOPX", "TSC22D1", "FAM30A", "MEF2C", "TPSB2", "GNG11", "TPT1", "TPSD1", "NRIP1", "ZBTB20", "KMT2A", "ST3GAL1", "RBPMS", "H1F0", "EMP1", "MEIS1")

Progenitor = c("CDK6", "HSP90AB1", "EEF1B2", "PCNP", "HINT1", "LRRC75A-AS1", "PEBP1", "LOC107984974", "H2AFY", "EEF1A1", "SMIM24", "PSME1", "SOX4", "EEF1G", "EBPL", "EIF4B", "PARP1", "MEST", "TFDP2", "ATP5G2", "NAP1L1", "TPM4", "SPN", "SPINK2", "SELL")

GMP = c("MPO", "SERPINB1", "CLEC11A", "CALR", "NUCB2", "CSF3R", "AZU1", "ELANE", "CD38", "FABP5", "LOC100419170", "HNRNPDL", "HSPB1", "RNA5-8S", "C12orf57", "TRH", "RUNX1T1", "POU4F1", "CEBPE", "LINC01835", "SNHG5", "THSD7A", "PRTN3")

Promono = c("SRGN", "CSTA", "MRPL33", "ANKRD28", "RNASE2", "PLD3", "RETN", "HSPA5", "EMB", "ZFR", "FAM107B", "MS4A3", "CTSG", "SLC44A1", "FUT4", "RNASE3")

Monocyte = c("S100A9", "S100A8", "S100A12", "PSAP", "FCN1", "SERPINA1", "CTSS", "DUSP1", "CEBPB", "NFKBIA", "VCAN", "CD14", "NAMPT", "TNFAIP2", "SLC11A1", "LILRB3")

cDC =  c("HLA-DQA1", "HLA-DRA", "HLA-DRB6", "HLA-DPB1", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "CST3", "DBI", "CRIP1", "CAP1", "TMSB10", "SAMHD1", "JAML", "ACTB", "NAPSB",  "CLEC10A", "GPX1", "ITGB7",  "FTH1P3", "FCER1A",  "S100B", "CPVL", "ALDH2", "CLEC4A")



all = c( HSC, Progenitor, GMP, Promono, Monocyte, cDC )
cell_cluster = rep(c("HSC-like","Progenitor-like","GMP-like","Promono-like","Monocyte-like","cDC-like"), times = c(length(HSC), length(Progenitor), length(GMP), length(Promono), length(Monocyte), length(cDC)))
names( cell_cluster ) = all
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds"

seurat_object = readRDS( rds_file )

# express_matrix_10000_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/backspin/express_matrix_10000.txt"
# express_matrix_10000 = read.table(  express_matrix_10000_file,header = TRUE )
# express_matrix_20000_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/backspin/express_matrix_20000.txt"
# express_matrix_20000 = read.table(  express_matrix_20000_file,header = TRUE )

# all_10000 = all[ all %in% rownames( express_matrix_10000 )  ]
# all_20000 = all[ all %in% rownames( express_matrix_20000 )  ]

# all_10000_express_matrix = express_matrix_10000[all_10000,]
# all_20000_express_matrix = express_matrix_20000[all_20000,]
# all_express_matrix = rbind( all_10000_express_matrix, all_20000_express_matrix )
# rowName = c( all_10000, all_20000 )
# cell_cluster = cell_cluster[rowName]
# cell_cluster = sort(cell_cluster)
# all_express_matrix = all_express_matrix[names(cell_cluster),]
# all_express_matrix_log = log(all_express_matrix+1)
# for (i in 1:dim(all_express_matrix_log)[1]){
# 	all_express_matrix_log[i,] = all_express_matrix_log[i,]/max(all_express_matrix_log[i,])
# }

express_matrix = seurat_object[["RNA"]]
all = all[ all %in% rownames( express_matrix ) ]
cell_cluster = cell_cluster[all]
#cell_cluster = sort(cell_cluster)
all_express_matrix = express_matrix[names(cell_cluster),]
all_express_matrix_log = all_express_matrix #数据已经为标准化之后，所以不需要进行log转换。

Sample = seurat_object[["sample"]]
FAB_class = seurat_object[["FAB_class"]]
response = seurat_object[["response"]]
cluster = as.character(seurat_object@active.ident)
names(cluster) = names(seurat_object@active.ident)


annotation_row = data.frame(CellCluster = cell_cluster)

fusion_gene=Sample$sample
fusion_gene[which(fusion_gene == "P101_M2_CR")] = "AML-ETO"
fusion_gene[which(fusion_gene == "P105_unknow_PR")] = "EV11"
fusion_gene[which(fusion_gene == "P106_M4_PR")] = "CBFB/MYH11"
fusion_gene[which(fusion_gene == "P107_M7_PRR")] ="EV11(+)"
fusion_gene[which(fusion_gene == "P108_M2_CR")] = "negative"
fusion_gene[which(fusion_gene == "P109_M3_unknow")] = "unknow"
fusion_gene[which(fusion_gene == "P110_M5_unknow")] = "FILT3-ITD"
fusion_gene[which(fusion_gene == "P111_M5_unknow")] = "negative"
fusion_gene[which(fusion_gene == "P112_unknow_unknow")] = "MLL/AF10"
names(fusion_gene) = rownames(Sample)

temp=Sample$sample
names(temp) = rownames(Sample)
temp = temp[cells]  #取一部分细胞
temp = sort(temp)
annotation_col = data.frame( Sample = temp )
all_express_matrix_log_plot = all_express_matrix_log[,names(temp)]
all_express_matrix_log_plot = t(scale(t(as.matrix(all_express_matrix_log_plot)),scale=FALSE))
plots = pheatmap( all_express_matrix_log_plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_row = annotation_row, annotation_col = annotation_col, color = colorRampPalette(colors = c("blue","grey","red"))(length(c(seq(-4,4,by=0.01)))),breaks= c(seq(-4,4,by=0.01)) )
ggsave(plots,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_AML/cell_AML_sample_all.png", height = 30,width = 20, limitsize = FALSE)

cluster = sort(cluster)
annotation_col = data.frame( Cluster = cluster )
all_express_matrix_log_plot = all_express_matrix_log[,names(cluster)]
plots = pheatmap( all_express_matrix_log_plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_row = annotation_row, annotation_col = annotation_col  )
ggsave(plots,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_AML/cell_AML_Cluster.png", height = 30,width = 20, limitsize = FALSE)

temp=FAB_class$FAB_class
names(temp) = rownames(FAB_class)
temp = temp[all_cells] 
temp = sort(temp)
fusion_gene = fusion_gene[names(temp)]
annotation_col = data.frame( FAB_class = temp, Fusion_gene = fusion_gene )
all_express_matrix_log_plot = all_express_matrix_log[,names(temp)]
all_express_matrix_log_plot =  t(scale(t(as.matrix(all_express_matrix_log_plot)),scale=FALSE))
plots = pheatmap( all_express_matrix_log_plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_row = annotation_row, annotation_col = annotation_col, color = colorRampPalette(colors = c("blue","grey","red"))(length(c(seq(-4,4,by=0.01)))),breaks= c(seq(-4,4,by=0.01))   )
ggsave(plots,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_AML/cell_AML_Cluster_FAB_class.png", height = 30,width = 20, limitsize = FALSE)


temp=response$response
names(temp) = rownames(response)
temp = sort(temp)
annotation_col = data.frame( response = temp )
all_express_matrix_log_plot = all_express_matrix_log[,names(temp)]
plots = pheatmap( all_express_matrix_log_plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_row = annotation_row, annotation_col = annotation_col  )
ggsave(plots,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_AML/cell_AML_Cluster_response.png", height = 30,width = 20, limitsize = FALSE)


##############提取肿瘤细胞绘制
all_cells = c()
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds"
seurat_object = readRDS( rds_file )
blast<-c( "LSC","HPC", "HSC" )
all_cluster = as.character(seurat_object@active.ident)
names(all_cluster) = colnames( seurat_object )
all_sample = seurat_object@meta.data$stim
names(all_sample) = colnames( seurat_object )
for (i in c( names(all_sample) )){
	if ( all_cluster[i] %in% blast){
		cells = c(cells,i)
		all_cells =c(all_cells,i)
	}
}


cells = c()
cellsType = c()
P101_cells = c()
P101_cellsType = c()
P101 <- readRDS("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P101SR_AllSample_CellCyleRegression_preprocess_removeDoublet/P101SR_26SampleAlignmentSeurat.rds")
blast<-c("stemCell_1","stemCell_2","stemCell_3","prog","GMP")
P101_cluster = as.character(P101@meta.data$orig.ident)
names(P101_cluster) = colnames( P101 )
P101_sample = P101@meta.data$stim
names(P101_sample) = colnames( P101 )
for (i in c( names(P101_sample) )){
	if ( P101_sample[i] == "P101SR00-PB" && P101_cluster[i] %in% blast){
		P101_cells =c(P101_cells,i)
		P101_cellsType = c(P101_cellsType, P101_cluster[i])
	}
}
names(P101_cellsType) = P101_cells
P101_cellsType = sort( P101_cellsType )
P101_cells = names( P101_cellsType )
cells = c( cells, P101_cells )
cellsType = c( cellsType, P101_cellsType )

P105_cells = c()
P105_cellsType = c()
P105 <-  readRDS("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P105SR_AllSample_CellCyleRegression_preprocess_removeDoublet/P105SR_38SampleAlignmentSeurat.rds")
blast<-c("prog_1","prog_2","prog_3","prog_4","cDC","CD14mono","stemCell","granu")
P105_cluster = as.character(P105@meta.data$orig.ident)
names(P105_cluster) = colnames( P105 )
P105_sample = P105@meta.data$stim
names(P105_sample) = colnames( P105 )
for (i in c( names(P105_sample) )){
	if ( P105_sample[i] == "P105SR00-PB" && P105_cluster[i] %in% blast){
		P105_cells =c(P105_cells,i)
		P105_cellsType = c(P105_cellsType, P105_cluster[i])
	}
	if ( P105_sample[i] == "P105SR00-BM" && P105_cluster[i] %in% blast){
		P105_cells =c(P105_cells,i)
		P105_cellsType = c(P105_cellsType, P105_cluster[i])
	}
}
names(P105_cellsType) = P105_cells
P105_cellsType = sort( P105_cellsType )
P105_cells = names( P105_cellsType )
cells = c( cells, P105_cells )
cellsType = c( cellsType, P105_cellsType )

P106_cells =c()
P106_cellsType = c()
P106 <-  readRDS("/asnas/wangqf_group/yuxx/pediatricAML/scRNASeqAnalysis/SeuratCellCycleRes/output/P106.combined_final_tutorial.rds")
blast<-c("HSC","CD38+", "IL3RA+", "unknown_1", "unknown_2")
P106_cluster = as.character(P106@active.ident)
P106_sample = as.character(P106@meta.data$orig.ident)
new_col_name = c()
raw_col_name = colnames( P106 )
for (i in 1:length(raw_col_name)){
	new_col_name = c(new_col_name,paste(P106_sample[i],strsplit(raw_col_name[i],"-")[[1]][1],sep = "_"))
}
names(P106_cluster) = new_col_name
names(P106_sample) = new_col_name
for (i in c( names(P106_sample) )){
	if ( P106_sample[i] == "P106SR00-PB" && P106_cluster[i] %in% blast){
		P106_cells =c(P106_cells,i)
		P106_cellsType = c(P106_cellsType, P106_cluster[i])
	}
	if ( P106_sample[i] == "P106SR00-BM" && P106_cluster[i] %in% blast){
		P106_cells =c(P106_cells,i)
		P106_cellsType = c(P106_cellsType, P106_cluster[i])
	}
}
names(P106_cellsType) = P106_cells
P106_cellsType = sort( P106_cellsType )
P106_cells = names( P106_cellsType )
cells = c( cells, P106_cells )
cellsType = c( cellsType, P106_cellsType )



########################################################cell_factor################################################
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
library(pheatmap)

interleukin = c("IL1B","IL1RN","IL2","IL2RA","IL3","IL4","IL6","IL10","IL12","IL12B","IL17A","IL18","IL23A","IL33","IL12A")
chemokine = c( "CCL3","CCL4","CCL5","CCL8","CCL7","CXCL1","CXCL12","CX3CL1","CXCL10","CXCL5","PF4","CXCL8")
tumor_necrosis_factor = c( "TNF","TNFRSF1A","TNFRSF1B","TNFSF10","FASLG","TNFRSF6B")
interferon = c("IFNA1","IFNA2","IFNB1","IFNG")
transforming_growth_factor = c( "TGFA","TGFB1","TGFB2","TGFB3","BMP2","BMP4","BMP6","BMP7")
growth_factor = c( "VEGFA","EGF","PDGFB","FGF7","FGF19","FGF10","FGF4","FGF12","FGF3","FGF6","HGF","IGF1","IGF2","LIF","NGF","OSM" )
colony_stimulating_factor = c("CSF3","CSF1","CSF2","KITLG","EPO")
others = c("SERPINE1", "TREM1", "TREM2", "PTX3", "CD40LG", "IL1RL1", "AGER", "IL6ST", "SELPLG", "CD44", "SELE", "SELP", "SELL", "GLP1R", "IRS1", "GCG", "LEP", "VSNL1", "BDNF", "NGF", "PLAT", "F3", "F9", "TNFRSF9", "CD27", "CD86", "CTLA4", "CD274", "PDCD1LG2", "PDCD1", "HAVCR2", "LAG3", "LGALS9")
adhesion_molecule = c("ICAM1", "ICAM2", "NCAM1", "ALCAM", "VCAM1", "PECAM1", "EPCAM", "ICAM3")

rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds"
seurat_object = readRDS( rds_file )

width = 6
height = 6

features = chemokine
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/chemokine.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = tumor_necrosis_factor
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/tumor_necrosis_factor.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = interferon
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/interferon.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = transforming_growth_factor
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/transforming_growth_factor.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = growth_factor
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/growth_factor.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = colony_stimulating_factor
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/colony_stimulating_factor.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = others
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/others3.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )

features = adhesion_molecule
print(length(features))
features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
print(length(features))

plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE, min.cutoff = "q9",split.by= "stim" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/cell_factor/adhesion_molecule.features.png", plots, height = height*length(features),limitsize = FALSE, width = 15*width  )
###################################################################################################################

###########################################绘制count数及gene数分布图################################################

express_matrix = read.table("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2017Science/expression_matrix_tpm.txt",header=TRUE,sep = "\t")

library(ggplot2)
paper_name = "2017Science"
gene_number = c()
for (i in 1:dim(express_matrix)[2]){
	gene_number = c(gene_number,sum(express_matrix[,i] >0))
}
gene_number = sort(gene_number)
pdf(paste("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/", paper_name, "/cell_number_histograms.pdf",sep=""))
plot(1:length(gene_number),gene_number,pch=20, xlab = "cell")
dev.off()
breaks = seq(0, 20000, 100)
breaks[length(breaks)] = max(max(gene_number),20000)
value = c(0,as.numeric(table(cut(gene_number,breaks = breaks))))
plots = plot( breaks, value, main = paper_name, xlab = "gene number", ylab = "cell number", pch = 20)
pdf(paste("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/", paper_name, "/cell_number.pdf",sep=""))
plot( breaks, value, main = paper_name, xlab = "gene number", ylab = "cell number", pch = 20)
dev.off()



cell_number = c()
for (i in 1:dim(express_matrix)[1]){
	cell_number = c(cell_number,sum(express_matrix[i,] >0))
}
cell_number = sort(cell_number)

pdf(paste("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/", paper_name, "/gene_number.pdf",sep=""))
plot( 1:length(cell_number), cell_numberbreaks, xlab = "gene", pch = 20)
dev.off()


breaks = seq(0, 1100, 50)
breaks[length(breaks)] = max(max(cell_number),1100)
value = c(0,as.numeric(table(cut(cell_number,breaks = breaks))))
plots = plot( breaks, value, main = paper_name, xlab = "cell number", ylab = "gene number", pch = 20)
pdf(paste("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/", paper_name, "/gene_number.pdf",sep=""))
plot( breaks, value, main = paper_name, xlab = "gene number", ylab = "cell number", pch = 20)
dev.off()


#########################################MRD_Marker##################################################
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
MRD_Marker = c( "NCAM1", "ANPEP", "IL3RA", "CD4", "FCGR1A", "PROM1", "CD7", "ITGA2B", "TFRC", "CSPG4" )
seurat_object = readRDS("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds")
MRD_Marker = MRD_Marker[MRD_Marker %in% rownames(seurat_object[["RNA"]]@data)]
DefaultAssay(object = seurat_object) <- "RNA"
sampleName = names(table(seurat_object@meta.data$orig.ident))
for (i in 1:length(sampleName)) {
    temp <- rownames(seurat_object@meta.data)[which(seurat_object$stim == sampleName[i])]
    single <- subset(seurat_object, cells = temp)
    p1 <- FeaturePlot(object = single, cells = temp, features = MRD_Marker, pt.size = 1.5, order=TRUE) #min.cutoff="q9",
    ggsave(paste("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/MRD_cluster_marker_",sampleName[i],".png",sep=""),p1, width = 40, height = 30)
    
}

######################################标记样本#########################################################
######3.0版本
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

#rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB_bk.rds"
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2019Cell/2.0/cell_Dim10_cluster.rds"

seurat_object = readRDS( rds_file )

samples = names(table(seurat_object[["stim"]][,1]))
for ( i in samples){
	result = rep( "NoCoverage", length( seurat_object$stim ) )
	names(result) = names(seurat_object$stim)
	result[names(seurat_object$stim[which(seurat_object$stim == i)])] = "WT"
	seurat_object[[paste( "sample_", i, sep = "" )]] <- result
}

p <- 1:(length(samples)+3)
p <- paste(rep("p",1),p,sep = "")
p1 <- DimPlot(object = seurat_object, reduction="umap",group.by = "stim", order=names(sort(table(seurat_object@meta.data$stim))), pt.size = 0.1) + NoLegend()
p2 <- DimPlot(object = seurat_object, reduction="umap",pt.size = 0.1, label = TRUE, label.size=1)  + NoLegend()
p3 <- DimPlot(object = seurat_object, reduction="umap",group.by = "Phase",pt.size = 0.1)  + NoLegend()
assign( p[1], p1 );assign( p[2], p2 );assign( p[3], p3 )
for ( i in 1:length(samples)){
	name = paste( "sample_", samples[i], sep = "" )
	order = c( "WT","NoCoverage" )
	cols = c( "#CCCCCC","#CC0033")
	temp = DimPlot(object = seurat_object, reduction="umap", group.by = name, order=order, cols=cols, pt.size = 0.1, cells=colnames(seurat_object) ) + labs( title = samples[i] ) + NoLegend()
	assign( p[i+3], temp )
}
text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
nrow = ceiling(length(p)/2)
eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2019Cell/2.0/sample_split.png", plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )

##################################
############################2.0版本
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R/OldVersion"))
library(ggplot2)
library(cowplot)
library(Matrix)

rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2019Cell/2.0/cell_Dim20_cluster.rds"

seurat_object = readRDS( rds_file )

samples = names(table(seurat_object@meta.data$stim))
for ( i in samples){
	result = rep( "NoCoverage", length( seurat_object@meta.data$stim ) )
	names(result) = seurat_object@cell.names
	result[which(seurat_object@meta.data$stim == i)] = "WT"
	add_data = data.frame( result )
	colnames(add_data) = paste( "sample_", i, sep = "" )
	seurat_object<-AddMetaData(object = seurat_object, metadata = add_data)
}




p <- 1:(length(samples))
p <- paste(rep("p",1),p,sep = "")
for ( i in 1:length(samples)){
	name = paste( "sample_", samples[i], sep = "" )
	order = c( "WT","NoCoverage" )
	cols = c( "#CCCCCC","#CC0033")
	temp = TSNEPlot(object =seurat_object, group.by=name,plot.order=c("WT","NoCoverage"),do.return = TRUE, pt.size = 0.1,colors.use=c("#CCCCCC","#CC0033")) + NoLegend()
	
	#temp = DimPlot(object = seurat_object, reduction="umap", group.by = name, order=order, cols=cols, pt.size = 0.1, cells=colnames(seurat_object) ) + labs( title = samples[i] ) + NoLegend()
	assign( p[i], temp )
}
text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
nrow = ceiling(length(p)/2)
eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2019Cell/2.0/sample_split.png", plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )

TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_PTPN11_check","Mut_PTPN11_check","N","NoCoverage_PTPN11_check"),do.return = TRUE, pt.size = 5,colors.use=c("#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"))

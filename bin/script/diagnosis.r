#点亮不同的样本的肿瘤样本并进行统计分布情况
library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
library(pheatmap)

tumor_cluster = c( "unknow", "Granulocytes", "classical Mono", "non-classical Mono", "cDC", "ProMono", "GMP", "HSC", "Progenitor"  )

###########################################################################################################################
#                                                                                                                         #
#                                                 每个样本的初诊肿瘤细胞ID                                                #
#                                                                                                                         #
###########################################################################################################################

#######基本信息#######
rds_file = c( "/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P105SR_AllSample_CellCyleRegression_preprocess_removeDoublet/P105SR_38SampleAlignmentSeurat.rds",
"/asnas/wangqf_group/yuxx/pediatricAML/scRNASeqAnalysis/SeuratCellCycleRes/output/P106.combined_all_final_tutorial.rds",
"/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P107SR_AllSample_CellCyleRegression_preprocess_removeDoublet/P107SR_34SampleAlignmentSeurat.rds"
"/asnas/wangqf_group/zhouzihan/scRNA/Seurat/NewIdent/P108SR_41SampleAlignmentSeurat.rds",
"/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P110SR_AllSample_CellCyleRegression_preprocess_removeDoublet/P110SR_43SampleAlignmentSeurat.rds",
"/asnas/wangqf_group/zhouzihan/scRNA/Seurat/P111/P111SR_45SampleAlignmentSeurat.rds",
"/asnas/wangqf_group/yuxx/pediatricAML/scRNASeqAnalysis/SeuratRes/P114/output/P114.combined_all_final_tutorial.rds",
"/asnas/wangqf_group/yuxx/pediatricAML/scRNASeqAnalysis/SeuratRes/P115/output/P115.combined_all_final_tutorial.rds",
"/asnas/wangqf_group/huangyue/seurat3/P116/P116SR_SampleAlignment_29cc.rds",
"/asnas/wangqf_group/yuxx/pediatricAML/scRNASeqAnalysis/SeuratRes/P117/output/P117.combined_all_final_tutorial.rds",
"/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/SingleSample_seurat3.0/P119SR00-BM/P119SR00-BM_26SampleAlignment.rds"
)
sample_name = c( "P105SR00-BM", "P106SR00-BM", "P107SR00-BM", "P108SR00-BM", "P110SR00-BM", "P111SR00-BM", "P114SR00-BM", "P115SR00-BM", "P116SR00-BM", "P117SR00-BM", "P119SR00-BM")

tumor_type = list(  c("prog_1","prog_2","prog_3","prog_4","cDC","CD14mono","stemCell","granu"),
					c("CD38+","HSC","IL3RA+","Gran","unknown_1","unknown_2"),
					c( "stemCell","CD14Mono","CD16Mono" ),
					c("stem_cell","prog"),
					c("stemCell","GMP","prog","CD14Mono","CD16Mono"),
					c("CD14+Mono","CD16+Mono","CD34+Unknown"),
					c("stem cells_1","stem cells_2","stem cells_3","unknown"),
					c("CMP","GMP-1","GMP-2","CD34-CD38+IL3RA+","unknown_2"),
					c("granu", "stem"),
					c("CMP","GMP","20","cDC","Classical Mono"),
					c("GMP","Granu","stemCell","Prog"))

###########输出每个样本的cellID#############

outdir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/tumor_cellID/"

for (i in 7:length( rds_file )){
	sample_rds_file = rds_file[i]
	sample = sample_name[i]
	blast = tumor_type[[i]]
	print(sample)
	cells = c()
	cellsType = c()
	seurat_object = readRDS( sample_rds_file )
	cluster = as.character(seurat_object@meta.data$orig.ident)
	names(cluster) = colnames( seurat_object )
	seurat_object_sample = seurat_object@meta.data$stim
	names(seurat_object_sample) = colnames( seurat_object )
	
	for (i in c( names(seurat_object_sample) )){
		if ( seurat_object_sample[i] == sample && cluster[i] %in% blast){
			cells =c(cells,i)
			cellsType = c(cellsType, cluster[i])
		}
	}
	print(length(cells))
	outfile = paste( outdir, sample, "_tumor_cellID.txt",sep = "" )
	write.table( cells, outfile, quote = FALSE, row.names = FALSE,  col.names = FALSE)
}

###########################################################################################################################
#                                                                                                                         #
#                                                      点亮细胞                                                           #
#                                                                                                                         #
###########################################################################################################################

##点亮样本的分布
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/diagnosis_Dim50_Res3.8_mark_mutation_and_celltype.rds"
seurat_object = readRDS( rds_file )

samples = names(table(seurat_object[["stim"]][,1]))
for ( i in 1:length(samples)){
	result = rep( "other sample", length( seurat_object$stim ) )
	names(result) = names(seurat_object$stim)
	result[names(seurat_object$stim[which(seurat_object$stim == samples[i])])] = samples[i]
	seurat_object[[paste( "sample_", samples[i], sep = "" )]] <- result
}

p <- 1:length(samples)
p <- paste(rep("p",1),p,sep = "")
col_plot = rainbow(length(samples)+1)
for ( i in 1:length(samples)){
	name = paste( "sample_", samples[i], sep = "" )
	order = c( samples[i] ,"other sample" )
	cols = c( "#CCCCCC", col_plot[i] )
	temp = DimPlot(object = seurat_object, reduction="umap", group.by = name, order=order, cols=cols, pt.size = 0.1, cells=colnames(seurat_object) ) + labs( title = samples[i] ) + NoLegend()
	assign( p[i], temp )
}
text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
nrow = ceiling(length(p)/2)
eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/sample.png", plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )


##点亮突变细胞
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/diagnosis_Dim50_Res3.8_mark_mutation_and_celltype.rds"

seurat_object = readRDS( rds_file )

sample_name = c( "P105SR00-BM", "P106SR00-BM", "P107SR00-BM", "P108SR00-BM", "P110SR00-BM", "P111SR00-BM", "P114SR00-BM", "P115SR00-BM", "P116SR00-BM", "P117SR00-BM", "P119SR00-BM")
tumor_cellID_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/tumor_cellID/"

tumor_cluster = c( "unknow", "Granulocytes", "classical Mono", "non-classical Mono", "cDC", "ProMono", "GMP", "HSC", "Progenitor"  )

cluster = as.character(seurat_object@active.ident)
names(cluster) = names(seurat_object@active.ident)

normal = rep( "other cells", length( seurat_object$stim ) )
names(normal) = names(seurat_object$stim)

for (i in 1:length(sample_name)){
	tumor_cellID_file = paste( tumor_cellID_dir, sample_name[i], "_tumor_cellID.txt", sep = "" )
	tumor_cellID = as.character( read.table( tumor_cellID_file )[,1] )
	tumor_cellID = names( cluster[tumor_cellID][which( cluster[tumor_cellID] %in% tumor_cluster)] )
	result = rep( "other cells", length( seurat_object$stim ) )
	names(result) = names(seurat_object$stim)
	result[ tumor_cellID ] = paste( sample_name[i], "malignant-like", sep = "_")
	normal[ tumor_cellID ] = "malignant-like"
	seurat_object[[paste( sample_name[i], "_tumor", sep = "" )]] <- result
}
seurat_object[["All_tumor"]] <- normal


p <- 1:(length(sample_name)+1)
p <- paste(rep("p",1),p,sep = "")
col_plot = rainbow(length(sample_name)+1)

for ( i in 1:length(sample_name)){
	name = paste( sample_name[i], "_tumor", sep = "" )
	order = c( paste( sample_name[i], "malignant-like", sep = "_") ,"other cells" )
	cols = c( "#CCCCCC", col_plot[i] )
	temp = DimPlot(object = seurat_object, reduction="umap", group.by = name, order=order, cols=cols, pt.size = 0.1, cells=colnames(seurat_object) ) + labs( title = sample_name[i] ) + NoLegend()
	assign( p[i], temp )
}

temp = DimPlot(object = seurat_object, reduction="umap", group.by = "All_tumor", order=c("malignant-like","other cells"), cols=c("#CCCCCC", col_plot[length(sample_name)+1]), pt.size = 0.1, cells=colnames(seurat_object) ) + labs( title = "All" ) + NoLegend()
assign( p[length(sample_name)+1], temp )

text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
nrow = ceiling(length(p)/2)
eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/sample_tumor.png", plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )


###########################################################################################################################
#                                                                                                                         #
#                                               统计每个样本中肿瘤细胞的比例                                              #
#                                                                                                                         #
###########################################################################################################################

rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/diagnosis_Dim50_Res3.8_mark_mutation_and_celltype.rds"
seurat_object = readRDS( rds_file )

sample_name = c( "P105SR00-BM", "P106SR00-BM", "P107SR00-BM", "P108SR00-BM", "P110SR00-BM", "P111SR00-BM", "P114SR00-BM", "P115SR00-BM", "P116SR00-BM", "P117SR00-BM", "P119SR00-BM")
tumor_cellID_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/tumor_cellID/"

all_cell_type = levels(seurat_object@active.ident)
result_count = matrix(0,nrow = length(sample_name), ncol = length(all_cell_type) )
rownames(result_count) = sample_name; colnames(result_count) = all_cell_type
result_rate = matrix(0,nrow = length(sample_name), ncol = length(all_cell_type) )
rownames(result_rate) = sample_name; colnames(result_rate) = all_cell_type

cluster = as.character(seurat_object@active.ident)
names(cluster) = names(seurat_object@active.ident)

for (i in 1:length(sample_name)){
	tumor_cellID_file = paste( tumor_cellID_dir, sample_name[i], "_tumor_cellID.txt", sep = "" )
	tumor_cellID = as.character( read.table( tumor_cellID_file )[,1] )
	tumor_cellID = names( cluster[tumor_cellID][which( cluster[tumor_cellID] %in% tumor_cluster)] )
	count = table( cluster[tumor_cellID] )
	result_count[sample_name[i], names(count) ] = count
	result_rate[ sample_name[i], ] = result_count[sample_name[i], ]/sum(result_count[sample_name[i], ])
}
write.table( result_count, "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/tumor_cell_count.txt", quote = FALSE, row.names = TRUE,  col.names = TRUE, sep= "\t")
result_rate<-round(result_rate,4)
write.table( result_rate, "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/tumor_cell_rate.txt", quote = FALSE, row.names = TRUE,  col.names = TRUE, sep= "\t")


###########################################################################################################################
#                                                                                                                         #
#                                                 cell 6个signature                                                       #
#                                                                                                                         #
###########################################################################################################################
HSC = c("HOPX", "TSC22D1", "FAM30A", "MEF2C", "TPSB2", "GNG11", "TPT1", "TPSD1", "NRIP1", "ZBTB20", "KMT2A", "ST3GAL1", "RBPMS", "H1F0", "EMP1", "MEIS1")

Progenitor = c("CDK6", "HSP90AB1", "EEF1B2", "PCNP", "HINT1", "LRRC75A-AS1", "PEBP1", "LOC107984974", "H2AFY", "EEF1A1", "SMIM24", "PSME1", "SOX4", "EEF1G", "EBPL", "EIF4B", "PARP1", "MEST", "TFDP2", "ATP5G2", "NAP1L1", "TPM4", "SPN", "SPINK2", "SELL")

GMP = c("MPO", "SERPINB1", "CLEC11A", "CALR", "NUCB2", "CSF3R", "AZU1", "ELANE", "CD38", "FABP5", "LOC100419170", "HNRNPDL", "HSPB1", "RNA5-8S", "C12orf57", "TRH", "RUNX1T1", "POU4F1", "CEBPE", "LINC01835", "SNHG5", "THSD7A", "PRTN3")

Promono = c("SRGN", "CSTA", "MRPL33", "ANKRD28", "RNASE2", "PLD3", "RETN", "HSPA5", "EMB", "ZFR", "FAM107B", "MS4A3", "CTSG", "SLC44A1", "FUT4", "RNASE3")

Monocyte = c("S100A9", "S100A8", "S100A12", "PSAP", "FCN1", "SERPINA1", "CTSS", "DUSP1", "CEBPB", "NFKBIA", "VCAN", "CD14", "NAMPT", "TNFAIP2", "SLC11A1", "LILRB3")

cDC =  c("HLA-DQA1", "HLA-DRA", "HLA-DRB6", "HLA-DPB1", "HLA-DRB1", "HLA-DQB1", "HLA-DPA1", "CST3", "DBI", "CRIP1", "CAP1", "TMSB10", "SAMHD1", "JAML", "ACTB", "NAPSB",  "CLEC10A", "GPX1", "ITGB7",  "FTH1P3", "FCER1A",  "S100B", "CPVL", "ALDH2", "CLEC4A")

all = c( HSC, Progenitor, GMP, Promono, Monocyte, cDC )
cell_cluster = rep(c("HSC-like","Progenitor-like","GMP-like","Promono-like","Monocyte-like","cDC-like"), times = c(length(HSC), length(Progenitor), length(GMP), length(Promono), length(Monocyte), length(cDC)))
names( cell_cluster ) = all

rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/diagnosis_Dim50_Res3.8_mark_mutation_and_celltype.rds"
seurat_object = readRDS( rds_file )

cluster = as.character(seurat_object@active.ident)
names(cluster) = names(seurat_object@active.ident)


sample_name = c( "P105SR00-BM", "P106SR00-BM", "P107SR00-BM", "P108SR00-BM", "P110SR00-BM", "P111SR00-BM", "P114SR00-BM", "P115SR00-BM", "P116SR00-BM", "P117SR00-BM", "P119SR00-BM")
sample_name = c( "P108SR00-BM", "P114SR00-BM", "P115SR00-BM", "P119SR00-BM", "P106SR00-BM","P116SR00-BM",  "P105SR00-BM", "P110SR00-BM", "P117SR00-BM", "P111SR00-BM", "P107SR00-BM")
tumor_cellID_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/tumor_cellID/"

cells = c()
sample = c()
cellType = c()
for (i in 1:length(sample_name)){
	tumor_cellID_file = paste( tumor_cellID_dir, sample_name[i], "_tumor_cellID.txt", sep = "" )
	tumor_cellID = as.character( read.table( tumor_cellID_file )[,1] )
	tumor_cellID = names( cluster[tumor_cellID][which( cluster[tumor_cellID] %in% tumor_cluster)] )
	tumor_cellID = names(sort(cluster[tumor_cellID]))
	cells = c(cells, tumor_cellID)
	sample = c( sample, rep( sample_name[i], length(tumor_cellID) ))
	cellType = c( cellType, cluster[tumor_cellID] )
}

express_matrix = seurat_object[["RNA"]]
all = all[ all %in% rownames( express_matrix ) ]
cell_cluster = cell_cluster[all]
all_express_matrix = express_matrix[names(cell_cluster),]
all_express_matrix_log = all_express_matrix[,cells]

annotation_col = data.frame( Sample = sample, cellType = cellType )
annotation_row = data.frame(CellCluster = cell_cluster)

all_express_matrix_log_plot = t(scale(t(as.matrix(all_express_matrix_log)),scale=FALSE))

plots = pheatmap( all_express_matrix_log_plot, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, annotation_row = annotation_row, annotation_col = annotation_col, color = colorRampPalette(colors = c("blue","grey","red"))(length(c(seq(-4,4,by=0.01)))),breaks= c(seq(-4,4,by=0.01)) )
ggsave(plots,file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/cell_6signature_scRNA.png", height = 30,width = 20, limitsize = FALSE)


###########################################################################################################################
#                                                                                                                         #
#                                                   寻找肿瘤群体中的各自的marker                                          #
#                                                                                                                         #
###########################################################################################################################

sample_name = c( "P105SR00-BM", "P106SR00-BM", "P107SR00-BM", "P108SR00-BM", "P110SR00-BM", "P111SR00-BM", "P114SR00-BM", "P115SR00-BM", "P116SR00-BM", "P117SR00-BM", "P119SR00-BM")
tumor_cellID_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/tumor_cellID/"

rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/diagnosis_Dim50_Res3.8_mark_mutation_and_celltype.rds"
seurat_object = readRDS( rds_file )

cluster = as.character(seurat_object@active.ident)
names(cluster) = names(seurat_object@active.ident)

result = cluster

cells = c()
cellType = c()
for (i in 1:length(sample_name)){
	tumor_cellID_file = paste( tumor_cellID_dir, sample_name[i], "_tumor_cellID.txt", sep = "" )
	tumor_cellID = as.character( read.table( tumor_cellID_file )[,1] )
	tumor_cellID = names( cluster[tumor_cellID][which( cluster[tumor_cellID] %in% tumor_cluster)] )
	cells = c(cells, tumor_cellID)
	cellType = c( cellType, cluster[tumor_cellID] )
}
names( cellType ) = cells
cellType_new = names(table( cellType ))
cellType = cellType[which(cellType %in% cellType_new)]
cellType_new = names(table(cellType))
for (i in cellType_new){
	result[names(cellType[which( cellType == i)])] = paste( i, "Like", sep = "" )
}

seurat_object[["Tumor_cell_type"]] <- result

###比较marker
Idents(object = seurat_object) = seurat_object@meta.data$Tumor_cell_type
outdir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/"

for (i in 1:length(cellType_new)){
	nk.markers = FindMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = paste( cellType_new[i], "Like", sep = "" ), ident.2 = cellType_new[i], grouping.by = "Tumor_cell_type")
	write.table(nk.markers ,file = paste( outdir,  paste( cellType_new[i], "Like", sep = "" ), "_", cellType_new[i],"_marker.txt", sep = ""), quote = FALSE,row.names = TRUE, sep = "\t")
}

tumor_cellType = paste( cellType_new, "Like", sep = "" )
for (i in 1:length(tumor_cellType)){
	nk.markers = FindMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, ident.1 = tumor_cellType[i], ident.2 = tumor_cellType[-i], grouping.by = "Tumor_cell_type")
	write.table(nk.markers ,file = paste( outdir,  tumor_cellType[i], "_other_tumor_marker.txt", sep = ""), quote = FALSE,row.names = TRUE, sep = "\t")
}

##############每个tumor cluster的gene list
marker_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/seurat_result_pc100/LightUp/"
tumor_cellType

all_gene_list = c()
marker_gene = list()
for (i in 1:length(tumor_cellType)){
	maker_file = paste( marker_dir, tumor_cellType[i], "_other_tumor_marker.txt", sep = "" )
	gene_matrix = read.table( maker_file, header = TRUE )
	temp = gene_matrix[ gene_matrix[,"p_val_adj"] < 0.01 , ]
	temp = temp[ temp[,"avg_logFC"] > 0 , ]
	marker_gene[[i]] = rownames(temp)
	all_gene_list = c(all_gene_list, rownames(temp))
	
}
names( marker_gene ) = tumor_cellType
uniqe_gene = names(table(all_gene_list)[table(all_gene_list) == 1])
for (i in tumor_cellType){
	print(i)
	print(length(marker_gene[[i]]))
	marker_gene[[i]] = marker_gene[[i]][which( marker_gene[[i]] %in% uniqe_gene )]
	marker_gene[[i]] = marker_gene[[i]][1:min(5,length(marker_gene[[i]]))]
	print(length(marker_gene[[i]]))
}

###########绘制TARGET数据

##TARGET_DATA_ANALYSIS.r脚本1-123行

tumor_cellType = c( "HSCLike", "GMPLike", "ProgenitorLike", "ProMonoLike", "classical MonoLike", "non-classical MonoLike", "cDCLike", "GranulocytesLike", "unknowLike" )

all = c()
cell_cluster = c()
for (i in tumor_cellType){
	all = c(all, marker_gene[[i]])
	cell_cluster = c( cell_cluster, rep( i, length(marker_gene[[i]]) ) )
}

names( cell_cluster ) = all

#express_matrix_bk = express_matrix
#clinicalData_matrix_bk = clinicalData_matrix
express_matrix = express_matrix_bk
clinicalData_matrix = clinicalData_matrix_bk



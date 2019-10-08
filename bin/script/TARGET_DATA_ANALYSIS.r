library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(dplyr)
library(ggplot2)
library(cowplot)
library(pheatmap)

option_list = list(
	make_option( c("-e", "--express_matrix"), type="character", default=NULL, 
				help="[Required] 表达矩阵（FPKM）", metavar = "character"),
	make_option( c("-o", "--outfile"), type="character", default=NULL, 
				help="[Required] 结果热图", metavar="character"),
	make_option(c("-f", "--feature"), type="character", default=NULL, 
				help="[Required] feature file", metavar="character"),
	make_option(c("-c", "--clinicalData"), type="character", default=NULL, 
				help="[Optional] 临床信息数据", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt$i = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/HTSeq_FPKM.merge.txt"
express_matrix = read.table( opt$i, head = TRUE, row.names = 1 )

# opt$f = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/mart_export.txt"
# features = read.table( opt$f, head = FALSE  )
# features = features[,c("Gene_stable_ID_version","Gene_name")]
# features_dup<-features[!duplicated(features, fromLast=TRUE), ]
# features_dup = as.matrix(features_dup)

opt$f = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/features.tsv"
features = read.table( opt$f, head = FALSE  )
features = features[,c(1,2)]
features_dup<-features[!duplicated(features, fromLast=TRUE), ]
features_dup = as.matrix(features_dup)

rowName = c()
raw_rowName = c()
for (i in rownames(express_matrix)){
	temp = strsplit( i, "\\.")[[1]][1]
	if ( temp %in% features_dup[,1 ]){
		rowName = c(rowName, features_dup[ features_dup[ ,1 ] == temp, 2 ])
		raw_rowName = c( raw_rowName,i )
	}
}

dup_names = names(table(rowName)[ table(rowName) > 1 ])
uniqueRowname = c()
unique_raw_Rowname = c()

for ( i in 1:length( rowName )){
	if ( !rowName[i] %in% dup_names){
		uniqueRowname = c( uniqueRowname,rowName[i] )
		unique_raw_Rowname = c( unique_raw_Rowname, raw_rowName[i] )
	}
}

all_express_matrix = express_matrix
express_matrix = express_matrix[ unique_raw_Rowname, ]
rownames( express_matrix ) = uniqueRowname

colName = colnames( express_matrix )
for ( i in 1:length(colName)){
	temp = strsplit( colName[i], "_" )[[1]]
	colName[i] = paste( temp[1:(length(temp)-1)], collapse = "_")
}
colnames( express_matrix ) = colName


opt$c = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/clinicalData.txt"
clinicalData = as.matrix( read.table( opt$c, head = TRUE, row.names = 1 ) )

all_clinicalData = clinicalData
clinicalData = clinicalData[colnames(express_matrix),]
clinicalData = t(clinicalData)

clinicalData_matrix = rep( "A", dim(clinicalData)[2] )
FAB = names( table(clinicalData[ "FAB", ]) )
rowName = c()
for ( i in FAB ){
	temp = rep( 0, dim(clinicalData)[2] )
	temp[ which(clinicalData[ "FAB", ] == i) ] = 1
	if ( i == "Unknown" ){
		i = "FAB_Unknown"
	}
	rowName = c( rowName, i )
	clinicalData_matrix = rbind( clinicalData_matrix, temp )
}

Fusion = names( table(clinicalData[ "Fusion", ]) )
for ( i in Fusion ){
	temp = rep( 0, dim(clinicalData)[2] )
	temp[ which(clinicalData[ "Fusion", ] == i) ] = 1
	if ( i == "Unknown" ){
		i = "Fusion_Unknown"
	}
	rowName = c( rowName, i )
	clinicalData_matrix = rbind( clinicalData_matrix, temp )
}

for (i in c("FLT3_ITD", "NPM", "CEBPA", "WT1", "c_Kit_Exon8", "c_Kit_Exon17")){
	temp = rep( 0, dim(clinicalData)[2] )
	temp[ which(clinicalData[ i, ] == "Yes") ] = 1
	temp[ which(clinicalData[ i, ] == "Unknown") ] = -1
	clinicalData_matrix = rbind( clinicalData_matrix, temp )
	rowName = c( rowName, i )
}

clinicalData_matrix = clinicalData_matrix[2:dim(clinicalData_matrix)[1],]
rownames( clinicalData_matrix ) = rowName
colnames( clinicalData_matrix ) = colnames( express_matrix )

################################根据fusion调整样本顺序###################################################################
Fusion = c( "CBFB_MYH11", "RUNX1_RUNX1T1", "RUNX1_CBFA2T2", "CBFA2T3_GLIS2", "KMT2A_AFF1", "KMT2A_ELL", "KMT2A_MLLT1","KMT2A_MLLT3", "KMT2A_MLLT4", "KMT2A_MLLT10", "KMT2A_SEPT6", "MLLT10_PICALM", "NUP98_HMGB3", "NUP98_HOXD13", "NUP98_KDM5A", "NUP98_NSD1", "NUP98_PHF23", "DEK_NUP214", "NPM1_MLF1", "MNX1_ETV6", "FUS_ERG", "FUS_FEV", "FUS_FLI1", "ETV6_INO8D", "negative" )

colName = c()
for (i in Fusion){
	colName = c( colName, colnames( clinicalData_matrix )[ which( clinicalData[ "Fusion", ] == i ) ] )
}

rowName = rownames(clinicalData_matrix)
rowName[10:34] = Fusion

clinicalData_matrix = clinicalData_matrix[rowName, ]

#########################################################################################################################

HSC = c("HOPX", "TSC22D1", "FAM30A", "MEF2C", "TPSB2", "GNG11", "TPT1", "TPSD1", "DST", "NRIP1", "ABCB1", "ZBTB20", "KMT2A", "ST3GAL1", "RBPMS", "TMEM74", "H1F0", "EMP1", "MEIS1", "PDZRN4", "XIRP2", "TMEM25", "C20orf203", "SLC6A13", "NPTX2", "ABCA9", "GABRA4", "CALCRL", "CLNK", "CRHBP")

Progenitor = c("CDK6", "HSP90AB1", "EEF1B2", "PCNP", "HINT1", "LRRC75A-AS1", "PEBP1", "LOC107984974", "H2AFY", "EEF1A1", "SMIM24", "PSME1", "SOX4", "EEF1G", "EBPL", "EIF4B", "PARP1", "MEST", "TMEM70", "TFDP2", "ATP5G2", "NAP1L1", "MSI2", "TPM4", "SPN", "SPINK2", "SELL", "TAPT1-AS1", "DSE", "LINC01623")

GMP = c("MPO", "SERPINB1", "CLEC11A", "AZU1", "CALR", "NUCB2", "CSF3R", "ELANE", "TRH", "RUNX1T1", "CD38", "POU4F1", "CEBPE", "LINC01835", "SNHG5", "FABP5", "LOC100419170", "HNRNPDL", "HSPB1", "RNA5-8S", "C12orf57", "THSD7A", "PRTN3", "CLEC5A", "PLPPR3", "IGFBP2", "PRRT4", "FBN2", "TSPOAP1", "FGFR1")

Promono = c("SRGN", "CSTA", "MRPL33", "ANKRD28", "RNASE2", "PLD3", "RETN", "HSPA5", "EMB", "ZFR", "FAM107B", "MS4A3", "CTSG", "SLC44A1", "SLPI", "FUT4", "RNASE3", "CCL23", "ATP8B4", "CLU", "KBTBD11", "PIWIL4", "DEFB1", "SERPINB10", "SESN3", "CD70", "PRLR", "LPL", "TP53INP2", "RNVU1-6")

Monocyte = c("S100A9", "S100A8", "PSAP", "FCN1", "S100A12", "SERPINA1", "BCL2A1", "FPR1", "CTSS", "DUSP1", "CEBPB", "NFKBIA", "VCAN", "CD14", "NAMPT", "TNFAIP2", "SLC11A1", "LILRB3", "BCL6", "CYP1B1", "MS4A10", "AQP9", "TLR4", "MAFB", "PLBD1", "THBS1", "C5AR1", "VNN2", "CR1", "APOBEC3A")

cDC =  c("HLA-DQA1", "HLA-DRA", "HLA-DRB6", "HLA-DPB1", "HLA-DRB1", "DBI", "CRIP1", "HLA-DQB1", "HLA-DPA1", "CAP1", "TMSB10", "SAMHD1", "JAML", "CST3", "ACTB", "NAPSB",  "CLEC10A", "GPX1", "ITGB7",  "FTH1P3", "FCER1A",  "MRC1", "HLA-DRB5", "S100B", "CPVL", "ALDH2", "CLEC4A", "PKIB", "HLA-DQA2", "HLA-DQB2")

all = c( HSC, Progenitor, GMP, Promono, Monocyte, cDC )
cell_cluster = rep(c("HSC-like","Progenitor-like","GMP-like","Promono-like","Monocyte-like","cDC-like"), times = c(length(HSC), length(Progenitor), length(GMP), length(Promono), length(Monocyte), length(cDC)))
names( cell_cluster ) = all

gene = intersect( all, rownames( express_matrix ) )

annotation_row = data.frame(CellCluster = cell_cluster[gene] )

plot_data = express_matrix[ gene, ]
#plot_data =  t(scale(t(as.matrix(plot_data)),scale=TRUE))
plot_data = log(plot_data+0.01)
plot_data =  t(scale(t(as.matrix(plot_data)),scale=FALSE))
plot_data = plot_data[ ,colName]


plots = pheatmap( plot_data, cluster_rows = FALSE, cluster_cols = TRUE, show_colnames = FALSE, annotation_row = annotation_row, color = colorRampPalette(colors = c("blue","grey","red"))(length(c(seq(-4,4,by=0.01)))),breaks= c(seq(-4,4,by=0.01)) )

opt$o = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/heatmap.png"

ggsave(plots,file = opt$o, height = 30,width = 20, limitsize = FALSE)


#sample_order = plots$tree_col$order
#clinicalData_matrix = clinicalData_matrix[,sample_order]
clinicalData_matrix = clinicalData_matrix[,colName]
clinicalData_matrix_plotdata=apply(clinicalData_matrix,2,as.numeric)
rownames(clinicalData_matrix_plotdata) = rownames( clinicalData_matrix )
plots = pheatmap( clinicalData_matrix_plotdata, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, color = colorRampPalette(colors = c("blue","white","grey"))(length(c(seq(-1,1,by=0.01)))),breaks= c(seq(-1,1,by=0.01)) )

opt$o = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/TARGET_analysis/rawData/heatmap_label.png"
ggsave(plots,file = opt$o, height = 10,width = 20, limitsize = FALSE)







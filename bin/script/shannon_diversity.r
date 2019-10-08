#用于计算单细胞数据中每个cluster的丰度
library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
usage <- function(){ #使用说明
	cat("\n使用说明\n")
	print_help(opt_parser)
}

option_list = list(
	make_option( c("-r", "--rds"), type="character", default=NULL, help="[Required] rds文件", metavar = "character"),
	make_option( c("-i", "--ini"), type="character", default=NULL, help="[Required] 配置文件", metavar = "character"),
	
	make_option( c("-o", "--out"), type="character", default=NULL, help="[Required] ", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

library(permute,lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(vegan,lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))

if(is.null( opt$i )){
	usage()
	q()
}
rds_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_bk.rds"
#seurat_object = readRDS( opt$i )
seurat_object = readRDS( rds_file )
cluster_num = as.character(seurat_object@meta.data$seurat_clusters)
names(cluster_num) = colnames(seurat_object)

Sample = seurat_object[["stim"]]$stim
names(Sample) = colnames(seurat_object)

cell_num_matrix = matrix( rep(0,length(table(Sample))*length(table(cluster_num))), nrow = length(table(cluster_num)), ncol = length(table(Sample)) )
colnames(cell_num_matrix) = names( table(Sample) )
rownames(cell_num_matrix) = names( table(cluster_num) )
for (i in rownames(cell_num_matrix)){
	temp = table(Sample[names(cluster_num[cluster_num == i])])
	cell_num_matrix[i,names(temp)] = temp
}
for (i in 1:dim(cell_num_matrix)[2]){
	sum_col = sum(cell_num_matrix[,i])
	for (j in 1:dim(cell_num_matrix)[1]){
		cell_num_matrix[j,i] = cell_num_matrix[j,i]/sum_col
	}
}

Shannon.Wiener <- diversity(cell_num_matrix, index = "shannon")
Shannon.Wiener = sort(Shannon.Wiener,decreasing = TRUE)

Sample.ini <- read.ini(opt$i)
CellType.ini <- Sample.ini$CellType
CellTypes <- names(CellType.ini)
boxplot_data <- list()
known = c()
for (i in CellTypes){
	print(i)
	print(CellType.ini[[i]])
	boxplot_data[[i]] = Shannon.Wiener[ strsplit(CellType.ini[[i]], ",")[[1]] ]
	known = c(known, strsplit(CellType.ini[[i]], ",")[[1]] )
}
boxplot_data[["unknown"]] = Shannon.Wiener[setdiff(names(Shannon.Wiener),known)]

temp = list()
mean_value= c()
for (i in names(boxplot_data)){
	mean_value = c(mean_value,median(boxplot_data[[i]]))
}
names(mean_value) = names(boxplot_data)
mean_value = sort(mean_value, decreasing = TRUE)
for (i in names(mean_value)){
	temp[[i]] = boxplot_data[[i]]
}
#boxplot(temp)

cellType = names(boxplot_data)
colors = rep(0,length(Shannon.Wiener));names(colors) = names(Shannon.Wiener)
colorsToSelect = rainbow(length(cellType))
for (i in 1:length(cellType)){
	colors[ names(boxplot_data[[cellType[i]]]) ] = colorsToSelect[i]
}

pdf( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/ShannonIndexBarplot.pdf",width=15 )
barplot(Shannon.Wiener,ylab="Shannon index", xlab = "cluster",col=colors,ylim=c(0,max(Shannon.Wiener)+0.2))
legend("topright",legend=cellType,col=colorsToSelect,ncol=3,bty="n",pch=15)
dev.off()

pdf( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/ShannonIndexBoxplot.pdf" )
boxplot(temp,col=rainbow(length(cellType)),par(las = "2"),cex.axis=0.6,ylab="Shannon index",outline = FALSE)
for (i in 1:length(temp)){
	points( rep(i,length(temp[[i]])), temp[[i]], pch=16 )
}
dev.off()




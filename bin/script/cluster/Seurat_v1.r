#/pnas/wangqf_group/zhangchen/hanxueshuai/R-3.5.0/bin/R
library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
########################################   optparse传参   ##############################################
option_list = list(
	make_option( c("-i", "--ini"), type="character", default=NULL, 
				help="配置文件，参见示例./example.ini", metavar = "character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, 
				help="结果输出目录", metavar="character"),
	make_option(c("-r", "--rds"), type="character", default=NULL, 
				help="读取保存的数据路径(.rds)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$i) || is.null(opt$o) ){ 
	cat("\n配置文档格式：\n")
	cat("
[Sample]
sample1=Abspath of cellranger's result of sample1, which includes barcodes.tsv, features.tsv(genes.tsv), matrix.mtx
sample2=Abspath of cellranger's result of sample2, which includes barcodes.tsv, features.tsv(genes.tsv), matrix.mtx
[Parameter]
features=features to plot. for example: CEACAM1,CEACAM6,CEACAM3,MTOR
min.features=200
min.nFeature_RNA=200
max.nFeature_RNA=2500
max.percent.mt=5
")
	cat("\n示例文档：/pnas/wangqf_group/suyx/PMO/pipeline/scRNA_seq_10X/bin/script/example.ini\n\n")
	print_help(opt_parser)
	q()
}

library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library( ggplot2 )
library(cowplot)

########################################   读取配置文件   #######################################
Sampe.ini <- read.ini(opt$i)
sample_name = names(Sampe.ini$Sample)
sample_dir_list = c()
for ( i in 1:length(sample_name) ){  #获取样本名字及路径
	temp = sample_name[i]
	sample_dir_list = c( sample_dir_list, Sampe.ini$Sample[i][[1]] )
}
features = strsplit(Sampe.ini$Parameter$features, ",")[[1]]
if (is.null( Sampe.ini$Parameter$min.features)){  #设置默认鉴定cell时最小features参数
	min.features = 200
}else{
	min.features = Sampe.ini$Parameter$min.features
}
if (is.null( Sampe.ini$Parameter$min.nFeature_RNA)){ #设置数据过滤时最小nFeature_RNA参数
	min.nFeature_RNA = 200
}else{
	min.nFeature_RNA = Sampe.ini$Parameter$min.nFeature_RNA
}
if (is.null( Sampe.ini$Parameter$max.nFeature_RNA)){ #设置数据过滤时最大nFeature_RNA参数
	max.nFeature_RNA = 2500
}else{
	max.nFeature_RNA = Sampe.ini$Parameter$max.nFeature_RNA
}
if (is.null( Sampe.ini$Parameter$max.percent.mt )){ #设置数据过滤时最大percent.mt参数
	max.percent.mt = 15
}else{
	max.percent.mt = Sampe.ini$Parameter$max.percent.mt
}

########################################   函数   ######################################################
cluster_calculation <- function( seurat_object, outfile, sample_name ){ #计算每个样本在每个cluster中的细胞数
	summary = c()
	clusters = names(table( seurat_object$seurat_clusters ))
	for (i in 1:length(clusters)){
		cluster_sample_cell_no = table(seurat_object$orig.ident[names(seurat_object$seurat_clusters[seurat_object$seurat_clusters == clusters[i]])])[sample_name]
		summary = c( summary, as.vector(cluster_sample_cell_no) ) 
	}
	summary_matrix = matrix( summary, ncol = length(sample_name), byrow = TRUE )
	rownames( summary_matrix ) <- clusters
	colnames( summary_matrix ) <- sample_name
	write.table(summary_matrix ,file = outfile, quote = FALSE,row.names = TRUE,sep = "\t" )
	return( summary_matrix )
}

qc_metrics_plots <- function( seurat_object, outdir, sample ){ #绘制qc指标
	plots = VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave( paste( outdir, "/", sample, "_QC_metrics1.pdf" ,sep = "" ), plots) 
	plot1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	plots = CombinePlots(plots = list(plot1, plot2))
	ggsave( paste( outdir, "/", sample, "_QC_metrics2.pdf" ,sep = "" ), plots )
}

variable_features_plots <- function( seurat_object, outdir, sample ){ #绘制特征分布图
	top10 <- head(VariableFeatures(seurat_object), 10)
	plot1 <- VariableFeaturePlot(seurat_object)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	plots = CombinePlots(plots = list(plot1, plot2))
	ggsave( paste( outdir, "/", sample, "_variable_features.pdf" ,sep = "" ) ,plots )
}

visualize_PCA <- function( seurat_object, outdir, sample ){ #绘制PCA分析结果
	plots = VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
	ggsave( paste( outdir, "/", sample, "_PCA_VizDimLoadings.pdf", sep = ""), plots )
	plots = DimPlot(seurat_object, reduction = "pca")
	ggsave( paste( outdir, "/", sample, "_PCA_DimPlot.pdf", sep = ""), plots )
	plots = DimHeatmap(seurat_object, dims = 1:3, cells = 500, balanced = TRUE)
	ggsave( paste( outdir, "/", sample, "_PCA_DimHeatmap.pdf", sep = ""), plots )
}

features_plot <- function( seurat_object, seurat_object.markers, features, outdir, sample ){ #绘制features图
	plots = VlnPlot(seurat_object, features = features, ncol=3)
	ggsave( paste( outdir, "/", sample, "_biomarkers_VlnPlot.pdf" ,sep = "" ), plots, height = 5*floor(length(features)/3),limitsize = FALSE )
	plots = FeaturePlot(seurat_object, features = features,ncol=2  )
	ggsave( paste( outdir, "/", sample, "_biomarkers_FeaturePlot.pdf" ,sep = "" ), plots, height = 5*floor(length(features)/2),limitsize = FALSE  )
	top10 <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
	plots = DoHeatmap(seurat_object, features = top10$gene) + NoLegend()
	ggsave( paste( outdir, "/", sample, "_biomarkers_DoHeatmap.pdf" ,sep = "" ), plots )
}

split_features <- function( seurat_object, features ){  #10个features为一组，以免文件过大，打开较慢
	all_features = rownames( seurat_object )
	result = list()
	temp = c()
	for (i in 1:length(features)){
		if (features[i] %in% all_features){
			temp = c( temp, features[i] )
			if (length(temp) == 10){
				result = c(result,list(temp))
				temp = c()
			}
		}
		
	}
	if (length(temp) > 0){
		result = c(result,temp)
	}
	return( result )
}

cell_no_barplot <- function( summary_matrix, outdir, sample ){
	clusters = rep( rownames( summary_matrix ), length(summary_matrix[1,]))
	sample_names = rep( colnames( summary_matrix ), each = length(summary_matrix[,1]))
	value = unlist(as.list(summary_matrix))
	ylabele_pos = summary_matrix[,1]
	for (i in 2:length(summary_matrix[1,])){
		ylabele_pos = c(ylabele_pos, rowSums(summary_matrix[, 1:i]) )
	}
	plot_data = data.frame( "cluster" = clusters, "sample" = sample_names, "value" = value, "ylabele_pos" = ylabele_pos )
	p1 <- ggplot(plot_data,aes(x = cluster,y = value, fill = sample))+geom_bar(stat='identity', position='stack')+geom_text(aes(y=ylabele_pos, label=value), vjust=1.6, color="white", size=3.5)+scale_fill_brewer(palette="Paired")+theme_minimal()
	p2 <- ggplot(plot_data,aes(x = cluster,y = value, fill = sample))+geom_bar(stat='identity', position='fill')
	plots = plot_grid(p1, p2)
	ggsave( paste( outdir, "/", sample, "_cellNo_BarPlot.pdf" ,sep = "" ), plots, width = 10 )
}
########################################   分析   #####################################################
seurat_list = c()
for (i in 1:length(sample_name)){  #读取cellranger参数并创建Seurat_object
	Seurat_object.data = Read10X( data.dir = sample_dir_list[i] )
	Seurat_object = CreateSeuratObject( counts = Seurat_object.data, project = sample_name[i], min.cells = 3,  min.features = min.features, assay = "RNA")
	Seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = Seurat_object, pattern = "^MT-")
	seurat_list = c(seurat_list, Seurat_object)
}

if (length( sample_name ) == 1 && is.null(opt$r)){
	i = 1
	seurat_object = seurat_list[[1]]
	qc_metrics_plots( seurat_object, opt$o, sample_name[i] )

	seurat_object <- subset(seurat_object, subset = nFeature_RNA > min.nFeature_RNA & nFeature_RNA < max.nFeature_RNA & percent.mt < max.percent.mt )
	seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
	seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
	variable_features_plots( seurat_object, opt$o, sample_name[i]  )
	
	all.genes <- rownames(seurat_object)
	seurat_object <- ScaleData(seurat_object, features = all.genes)
	seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
	visualize_PCA( seurat_object, opt$o, sample_name[i]  )
	
	seurat_object <- JackStraw(seurat_object, num.replicate = 100)
	seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
	plots = JackStrawPlot(seurat_object, dims = 1:15)
	ggsave( paste( opt$o, "/", sample_name[i], "_JackStrawPlot.pdf" ,sep = "" ), plots )
	plots = ElbowPlot(seurat_object)
	ggsave( paste( opt$o, "/", sample_name[i], "_ElbowPlot.pdf" ,sep = "" ), plots )
	
	seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
	seurat_object <- FindClusters(seurat_object, resolution = 0.5)
	seurat_object <- RunUMAP(seurat_object, dims = 1:10)
	plots = DimPlot(seurat_object, reduction = "umap")
	ggsave( paste( opt$o, "/", sample_name[i], "_umap.pdf" ,sep = "" ), plots )
	saveRDS(seurat_object, file = paste( opt$o, "/", sample_name[i], ".rds" ,sep = ""))

	seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	features_to_plot <- split_features( seurat_object, features )
	for (i in 1:length( features_to_plot )){
		features_plot( seurat_object, seurat_object.markers, features_to_plot[[i]], opt$o, paste(sample_name[i],"_",i,sep="") )
	}
}else if(length( sample_name ) == 1 && !is.null(opt$r)){
	i = 1
	seurat_object = readRDS( file = opt$r )
	seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	features_to_plot <- split_features( seurat_object, features )
	for (i in 1:length( features_to_plot )){
		features_plot( seurat_object, seurat_object.markers, features_to_plot[[i]], opt$o, paste(sample_name[i],"_",i,sep="") )
	}
}else{
	file_name = paste(sample_name, collapse = "_")
	for (i in 1:length(seurat_list)){
		seurat_list[[i]]$stim <- sample_name[i]
		seurat_list[[i]] = subset(seurat_list[[i]], subset = nFeature_RNA > min.nFeature_RNA)
		seurat_list[[i]] = NormalizeData( seurat_list[[i]], verbose = FALSE)
		seurat_list[[i]] = FindVariableFeatures( seurat_list[[i]], selection.method = "vst", nfeatures = 2000 )
	}
	seurat_list.anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:20)
	seurat_list.combined <- IntegrateData(anchorset = seurat_list.anchors, dims = 1:20)
	
	DefaultAssay(seurat_list.combined) <- "integrated"
	# Run the standard workflow for visualization and clustering
	seurat_list.combined <- ScaleData(seurat_list.combined, verbose = FALSE)
	seurat_list.combined <- RunPCA(seurat_list.combined, npcs = 30, verbose = FALSE)
	# t-SNE and Clustering
	seurat_list.combined <- RunUMAP(seurat_list.combined, reduction = "pca", dims = 1:20)
	seurat_list.combined <- FindNeighbors(seurat_list.combined, reduction = "pca", dims = 1:20)
	seurat_list.combined <- FindClusters(seurat_list.combined, resolution = 0.5)
	# Visualization
	p1 <- DimPlot(seurat_list.combined, reduction = "umap", group.by = "stim")
	p2 <- DimPlot(seurat_list.combined, reduction = "umap", label = TRUE)
	plots = plot_grid(p1, p2)
	ggsave( paste( opt$o, "/", file_name, "_cluster_DimPlot.pdf" ,sep = "" ), plots, width = 10 )
	plots = DimPlot(seurat_list.combined, reduction = "umap", split.by = "stim")
	ggsave( paste( opt$o, "/", file_name, "_sample_DimPlot.pdf" ,sep = "" ), plots, width = 5*length(sample_name) )
	DefaultAssay(seurat_list.combined) <- "RNA"
	nk.markers = FindAllMarkers(seurat_list.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, grouping.var = "stim")
	write.table(nk.markers ,file = paste( opt$o, "/", file_name, "_all_marker.txt" ,sep = "" ),quote = FALSE,row.names = TRUE,sep = "\t")
	summary_matrix <- cluster_calculation( seurat_list.combined, paste( opt$o, "/", file_name, "_summary_matrix.txt" ,sep = "" ), sample_name)
	cell_no_barplot( summary_matrix, opt$o, file_name )
	features_to_plot <- split_features( seurat_list.combined, features )
	for (i in 1:length(features_to_plot)){
		plots = FeaturePlot(seurat_list.combined, features = features_to_plot[[i]], min.cutoff = "q9",ncol = 2)
		ggsave( paste( opt$o, "/", file_name, "_",i, "_FeaturePlot.pdf" ,sep = "" ), plots,height = 5*floor(length(features_to_plot[[i]])/2),limitsize = FALSE, width = 10  )
	}
	saveRDS(seurat_list.combined, file = paste( opt$o, "/", file_name, ".rds" ,sep = ""))
}


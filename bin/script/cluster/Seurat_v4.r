#/pnas/wangqf_group/zhangchen/hanxueshuai/R-3.5.0/bin/R
library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)

#=======================================      函数    ================================================#

usage <- function(){ #使用说明
	cat("\n配置文档格式：\n")
	cat("
[Sample]
sample1=Abspath of cellranger's result of sample1, which includes barcodes.tsv, features.tsv(genes.tsv), matrix.mtx
sample2=Abspath of cellranger's result of sample2, which includes barcodes.tsv, features.tsv(genes.tsv), matrix.mtx
[Parameter]
#dim=25,26,27,28
features=features to plot. for example: CEACAM1,CEACAM6,CEACAM3,MTOR
min.features=200
min.nFeature_RNA=200
max.nFeature_RNA=25000
max.percent.mt=15
各参数含义请见示例文档
")
	cat("\n示例文档：/pnas/wangqf_group/suyx/PMO/pipeline/scRNA_seq_10X/bin/script/example.ini\n\n")
	print_help(opt_parser)
}

load_cellranger_result <- function( sample_name, sample_dir_list, outdir, max.percent.mt, min.nFeature_RNA, max.nFeature_RNA, min.features, file_name, datatype ){ #读取cellranger分析结果
	seurat_list = c()
	infor2count<-c("rawCellNo","remainCellNo","filtered",paste("percent.mito>=",max.percent.mt,"%",sep=""), paste("nGene<=",min.nFeature_RNA,sep=""), paste("nGene>=", max.nFeature_RNA, sep="")  )
	count<-matrix(0,nrow=length(sample_name),ncol=length(infor2count))
	rownames(count)<-sample_name
	colnames(count)<-infor2count
	for (i in 1:length(sample_name)){  #读取cellranger参数并创建Seurat_object
		print(sample_dir_list[i])
		if (datatype == "bc_matrix"){
			Seurat_object.data = Read10X( data.dir = sample_dir_list[i] )
			Seurat_object = CreateSeuratObject( counts = Seurat_object.data, project = sample_name[i], min.cells = 3,  min.features = min.features)
		}else{
			temp = as.matrix(read.table( file = sample_dir_list[i], as.is=TRUE, row.names=1, header=TRUE ))
			Seurat_object = CreateSeuratObject( counts = temp )
		}
		print( dim(Seurat_object) )
		Seurat_object<- RenameCells(Seurat_object,add.cell.id=sample_name[i])
		Seurat_object$stim <- sample_name[i]
		Seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = Seurat_object, pattern = "^MT-")
		seurat_list = c(seurat_list, Seurat_object)
		rawCellNo=dim(Seurat_object)[2]
		remainCellNo=dim(subset(Seurat_object@meta.data,subset = nFeature_RNA > min.nFeature_RNA & nFeature_RNA < max.nFeature_RNA & percent.mt < max.percent.mt ))[1]
		filtered = rawCellNo - remainCellNo
		count[sample_name[i],] = c( rawCellNo, remainCellNo,filtered,
									dim(subset(Seurat_object@meta.data,percent.mt>=max.percent.mt))[1],
									dim(subset(Seurat_object@meta.data,nFeature_RNA<=min.nFeature_RNA))[1],
									dim(subset(Seurat_object@meta.data,nFeature_RNA>=max.nFeature_RNA))[1])
	}
	write.table(count,paste(outdir, "/", file_name,"_SampleFilteringCount.txt",sep=""),sep="\t",quote = FALSE)
	return( seurat_list )
}

qc_metrics_plots <- function( seurat_object, outdir, sample ){ #绘制qc指标
	plots = VlnPlot(object = seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave( paste( outdir, "/", sample, "_QC_metrics1.pdf" ,sep = "" ), plots, width = 10,limitsize = FALSE ) 
	plot1 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	plots = CombinePlots(plots = list(plot1, plot2))
	ggsave( paste( outdir, "/", sample, "_QC_metrics2.pdf" ,sep = "" ), plots, width = 10,limitsize = FALSE )
}

filter_seurat_object <- function( seurat_object, max.percent.mt, min.nFeature_RNA, max.nFeature_RNA  ){ #过滤nFeature_RNA percent.mt并标准化数据
	seurat_object <- subset(seurat_object, subset = nFeature_RNA > min.nFeature_RNA & nFeature_RNA < max.nFeature_RNA & percent.mt < max.percent.mt )
	seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
	seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
	return( seurat_object )
}

variable_features_plots <- function( seurat_object, outdir, sample ){ #绘制特征分布图
	top10 <- head(VariableFeatures(seurat_object), 10)
	plot1 <- VariableFeaturePlot(seurat_object)
	plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	plots = CombinePlots(plots = list(plot1, plot2))
	ggsave( paste( outdir, "/", sample, "_variable_features.png" ,sep = "" ) ,plots,limitsize = FALSE, height = 10, width = 20 )
}

visualize_PCA <- function( seurat_object, outdir, sample ){ #绘制PCA分析结果
	plots = DimHeatmap(seurat_object, dims = 1:25, cells = 1000, balanced = TRUE, ncol = 4)
	ggsave( paste( outdir, "/", sample, "_PCA_DimHeatmap.png", sep = ""), plots, width = 20, height = 5*7,limitsize = FALSE  )

	plots = ElbowPlot(seurat_object, ndims = 30)
	ggsave( paste( opt$o, "/", sample, "_ElbowPlot.pdf" ,sep = "" ), plots,limitsize = FALSE )
}

sample_cluster_DimPlot <- function( seurat_object, outdir, sample, res, sample_name ){
	width = 6
	height = 6
	p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "stim", pt.size = 0.1 )
	p2 <- DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.1 )
	plots = plot_grid(p1, p2)
	ggsave( paste( outdir, "/", sample, "_cluster_DimPlot.pdf" ,sep = "" ), plots, width = width, limitsize = FALSE , height = height)
	plots <- DimPlot(seurat_object, reduction = "umap", group.by = "stim", width = width, height = height ) + NoLegend()
	ggsave( paste( outdir, "/", sample, "_cluster_DimPlot_NoLegend.pdf" ,sep = "" ), plots, width = width, limitsize = FALSE, height = height )
	plots = DimPlot(seurat_object, reduction = "umap", split.by = "stim")
	ggsave( paste( outdir, "/", sample, "_sample_DimPlot.pdf" ,sep = "" ), plots, width = width*2, limitsize = FALSE, height = height*2 )
	plots = DimPlot(seurat_object, reduction = "umap", split.by = "stim", group.by = "stim", ncol = 3)
	ggsave( paste( outdir, "/", sample, "_sample_DimPlot_bySample.pdf" ,sep = "" ), plots, width = width*2, limitsize = FALSE, height = height*2 )
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
}

features_plot <- function( seurat_object, seurat_object.markers, features, outdir, sample ){ #绘制features图
	width = 5
	height = 5
	plots = VlnPlot(seurat_object, features = features, ncol=4, pt.size = 0.1)
	ggsave( paste( outdir, "/", sample, "_biomarkers_VlnPlot.png" ,sep = "" ), plots, height = height*ceiling(length(features)/4),limitsize = FALSE, width = 4*width   )
	plots = FeaturePlot(seurat_object, features = features,ncol=4, pt.size = 0.1, order = TRUE )
	ggsave( paste( outdir, "/", sample, "_biomarkers_FeaturePlot.png" ,sep = "" ), plots, height = height*ceiling(length(features)/4),limitsize = FALSE, width = 4*width  )
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

Cluster_cellType <- function( CellType.ini, cluster ){
	CellTypes <- names(CellType.ini)
	new.cluster.ids <- rep(0,length( cluster ))
	all_CellType_cluster <- c()
	for ( i in 1:length(CellTypes)){
		CellType_name = CellType.ini[i]
		CellType_cluster = strsplit(CellType.ini[i][[1]], "," )[[1]]
		new.cluster.ids[ which( cluster %in% CellType_cluster )] <- CellTypes[i]
		all_CellType_cluster = c( all_CellType_cluster, CellType_cluster)
	}
	if (!all( all_CellType_cluster %in%  cluster)){
		temp = paste(all_CellType_cluster[ which( !(all_CellType_cluster %in%  cluster) ) ], collapse = ",")
		print( paste( "cluster ", temp, " 不在聚类结果中，请修改配置文件" ) )
		q()
	}
	if ( !all( table(all_CellType_cluster) == 1 ) ){
		temp = table(all_CellType_cluster)
		temp = paste(temp[which(temp > 1)], collapse = ",")
		print( paste( "cluster ", temp, " 出现在多个celltype中，请修改配置文件" ) )
		q()
	}
	new.cluster.ids[ which(new.cluster.ids ==  0)] = "unknow"
	return( new.cluster.ids )
}

cluster_calculation <- function( seurat_object, outfile, sample_name ){ #计算每个样本在每个cluster中的细胞数
	summary = c()
	clusters = names(table( seurat_object@active.ident ))
	for (i in 1:length(clusters)){
		cluster_sample_cell_no = table(seurat_object$orig.ident[names(seurat_object@active.ident[seurat_object@active.ident == clusters[i]])])[sample_name]
		summary = c( summary, as.vector(cluster_sample_cell_no) ) 
	}
	summary_matrix = matrix( summary, ncol = length(sample_name), byrow = TRUE )
	rownames( summary_matrix ) <- clusters
	colnames( summary_matrix ) <- sample_name
	write.table(summary_matrix ,file = outfile, quote = FALSE,row.names = TRUE,sep = "\t" )
	return( summary_matrix )
}

cell_no_barplot <- function( summary_matrix, outdir, sample ){ #绘制每个细胞在每个cluster中的细胞数
	summary_matrix = round(summary_matrix,2)
	clusters = rep( rownames( summary_matrix ), length(summary_matrix[1,]))
	sample_names = rep( colnames( summary_matrix ), each = length(summary_matrix[,1]))
	value = unlist(as.list(summary_matrix))
	ylabele_pos = c()
	for (i in 1:(length(summary_matrix[1,])-1)){
		ylabele_pos = c(ylabele_pos, rowSums(summary_matrix[, i:length(summary_matrix[1,])] ))
	}
	ylabele_pos = c(ylabele_pos, summary_matrix[, length(summary_matrix[1,])])
	plot_data = data.frame( "cluster" = clusters, "sample" = sample_names, "value" = value, "ylabele_pos" = ylabele_pos )
	p1 <- ggplot(plot_data,aes(x = cluster,y = value, fill = sample))+geom_bar(stat='identity', position='stack')+geom_text(aes(y=ylabele_pos, label=value), vjust=1.6, color="white", size=3.5)+scale_fill_brewer(palette="Paired")+theme_minimal()
	p2 <- ggplot(plot_data,aes(x = cluster,y = value, fill = sample))+geom_bar(stat='identity', position='fill')
	plots = plot_grid( p1, p2 )
	ggsave( paste( outdir, "/", sample, "_cellNo_BarPlot.pdf" ,sep = "" ), plots, width = 10,limitsize = FALSE )
}

mutation_marker_plot <- function( Mutation, seurat_object, outfile, res ){ #绘制突变图
	mutation_sample_gene <- names( Mutation )
	p <- 1:(length(mutation_sample_gene)+3)
	p <- paste(rep("p",length(res)),p,sep = "")
	p1 <- DimPlot(object = seurat_object, reduction="umap",group.by = "stim", order=names(sort(table(seurat_object@meta.data$stim))), pt.size = 0.1)
	p2 <- DimPlot(object = seurat_object, reduction="umap",pt.size = 0.1, label = TRUE, label.size=1)
	p3 <- DimPlot(object = seurat_object, reduction="umap",group.by = "Phase",pt.size = 0.1) 
	assign( p[1], p1 );assign( p[2], p2 );assign( p[3], p3 )
	
	for ( i in 1:length(mutation_sample_gene)){
		sample_gene = mutation_sample_gene[i]
		dir <- Mutation[[sample_gene]]
		gene_mutation = read.table( dir )
		sample = strsplit( sample_gene, "_" )[[1]][1]
		gene <- strsplit( sample_gene, "_" )[[1]][2]
		name = paste( "finalCellMutaionStatus", gene, sep = "_" )
		#cell.use = intersect(rownames(gene_mutation),rownames(seurat_object@meta.data[[name]]))
		
		seurat_object_sample = as.list(seurat_object[["stim"]])[[1]]
		names(seurat_object_sample) = rownames(seurat_object[["stim"]])
		cell.use = intersect(names(seurat_object_sample[seurat_object_sample == sample]),rownames(seurat_object@meta.data[[name]]))
		
		order = names( table( seurat_object[[name]][cell.use,] ) )
		if (length(order) == 4){
			order = c( "Mut","WT","further2confirm","NoCoverage" )
			cols = c( "#CCCCCC","#FFFF00","#FFCCCC","#CC0033")
		}else{
			order = c("Mut","WT","NoCoverage")
			cols=c("#CCCCCC","#FFFF00","#CC0033")
		}
		
		temp = DimPlot(object = seurat_object, reduction="umap",group.by = name, order=order, cols=cols,pt.size = 0.1, cells = cell.use) + labs( title = sample_gene )
		assign( p[i+3], temp )
	}
	text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
	nrow = ceiling(length(p)/2)
	eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))

	ggsave( outfile, plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )
}

mutation_marker_plot_NoLegend <- function( Mutation, seurat_object, outfile, res ){ #绘制突变图
	mutation_sample_gene <- names( Mutation )
	p <- 1:(length(mutation_sample_gene)+3)
	p <- paste(rep("p",length(res)),p,sep = "")
	p1 <- DimPlot(object = seurat_object, reduction="umap",group.by = "stim", order=names(sort(table(seurat_object@meta.data$stim))), pt.size = 0.1) + NoLegend()
	p2 <- DimPlot(object = seurat_object, reduction="umap",pt.size = 0.1, label = TRUE, label.size=1)  + NoLegend()
	p3 <- DimPlot(object = seurat_object, reduction="umap",group.by = "Phase",pt.size = 0.1)  + NoLegend()
	assign( p[1], p1 );assign( p[2], p2 );assign( p[3], p3 )
	
	for ( i in 1:length(mutation_sample_gene)){
		sample_gene = mutation_sample_gene[i]
		dir <- Mutation[[sample_gene]]
		gene_mutation = read.table( dir )
		sample = strsplit( sample_gene, "_" )[[1]][1]
		gene <- strsplit( sample_gene, "_" )[[1]][2]
		name = paste( "finalCellMutaionStatus", gene, sep = "_" )
		#cell.use = intersect(rownames(gene_mutation),rownames(seurat_object@meta.data[[name]]))
		
		seurat_object_sample = as.list(seurat_object[["stim"]])[[1]]
		names(seurat_object_sample) = rownames(seurat_object[["stim"]])
		cell.use = intersect(names(seurat_object_sample[seurat_object_sample == sample]),rownames(seurat_object@meta.data[[name]]))
		
		order = names( table( seurat_object[[name]][cell.use,] ) )
		if (length(order) == 4){
			order = c( "Mut","WT","further2confirm","NoCoverage" )
			cols = c( "#CCCCCC","#FFFF00","#FFCCCC","#CC0033")
		}else{
			order = c("Mut","WT","NoCoverage")
			cols=c("#CCCCCC","#FFFF00","#CC0033")
		}
		
		temp = DimPlot(object = seurat_object, reduction="umap",group.by = name, order=order, cols=cols, pt.size = 0.1, cells = cell.use) + labs( title = sample_gene ) + NoLegend()
		assign( p[i+3], temp )
	}
	text = paste( paste( "get(p[",1:length(p),"]) ", sep = "" ), collapse = ",")
	nrow = ceiling(length(p)/2)
	eval(parse(text = paste( "plots <- plot_grid( ", text, ", nrow=nrow,ncol=2 )",sep = "" )))

	ggsave( outfile, plots, height = 5*ceiling(length(p)/2),width=10,limitsize = FALSE )
}

#=======================================       分析       ================================================#

#=======================================    optparse传参  ================================================#
option_list = list(
	make_option( c("-i", "--ini"), type="character", default=NULL, 
				help="[Required] 配置文件，参见示例./example.ini", metavar = "character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, 
				help="[Required] 结果输出目录", metavar="character"),
	make_option(c("-d", "--datatype"), type="character", default="bc_matrix", 
				help="[Required] 数据类型，是表达矩阵还是cellranger分析结果。bc_matrix:cellranger分析结果 express_matrix:表达矩阵", metavar="character"),
	make_option(c("-r", "--rds"), type="character", default=NULL, 
				help="[Optional] 读取保存的数据路径(.rds)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (is.null(opt$i) || is.null(opt$o) ){ 
	usage()
	q()
}else if(!file.exists( opt$o )) {
	dir.create(opt$o)
}


#=======================================   读取配置文件   ================================================#
Sample.ini <- read.ini(opt$i)
sample_name = names(Sample.ini$Sample)
sample_dir_list = c()
for ( i in 1:length(sample_name) ){  #获取样本名字及路径
	temp = sample_name[i]
	sample_dir_list = c( sample_dir_list, Sample.ini$Sample[i][[1]] )
}

if ( is.null(Sample.ini$Analysis$analysis) || !( Sample.ini$Analysis$analysis %in% c("load_data_PCA", "cluster", "mark_mutation_and_celltype") ) ){ #确定分析类型是否正确
	print("请在配置文件中配置分析类型[Analysis]，参见示例文档")
	usage()
	q()
}

if (is.null( Sample.ini$Analysis$out_file )){
	file_name <- paste(sample_name, collapse = "_")
}else{
	file_name <- Sample.ini$Analysis$out_file
}

if ( Sample.ini$Analysis$analysis == "load_data_PCA" ){ #读取数据并分析至PCA部分
	if ( is.null(Sample.ini$Load_Data_PCA_Parameter$s.genes) || is.null(Sample.ini$Load_Data_PCA_Parameter$g2m.genes) ){ #读取配置文件中的参数
		print("请在配置文件中配置s.genes、g2m.genes参个参数，参见示例文档")
		usage()
		q()
	}else{
		s.genes = strsplit(Sample.ini$Load_Data_PCA_Parameter$s.genes, ",")[[1]]
		g2m.genes = strsplit(Sample.ini$Load_Data_PCA_Parameter$g2m.genes, ",")[[1]]
	}
	if (is.null( Sample.ini$Load_Data_PCA_Parameter$min.features)){  #设置默认鉴定cell时最小features参数
		min.features = 200
	}else{
		min.features = as.numeric( Sample.ini$Load_Data_PCA_Parameter$min.features )
	}
	if (is.null( Sample.ini$Load_Data_PCA_Parameter$min.nFeature_RNA)){ #设置数据过滤时最小nFeature_RNA参数
		min.nFeature_RNA = 200
	}else{
		min.nFeature_RNA = as.numeric(Sample.ini$Load_Data_PCA_Parameter$min.nFeature_RNA)
	}
	if (is.null( Sample.ini$Load_Data_PCA_Parameter$max.nFeature_RNA)){ #设置数据过滤时最大nFeature_RNA参数
		max.nFeature_RNA = 25000
	}else{
		max.nFeature_RNA = as.numeric(Sample.ini$Load_Data_PCA_Parameter$max.nFeature_RNA)
	}
	if (is.null( Sample.ini$Load_Data_PCA_Parameter$max.percent.mt )){ #设置数据过滤时最大percent.mt参数
		max.percent.mt = 15
	}else{
		max.percent.mt = as.numeric(Sample.ini$Load_Data_PCA_Parameter$max.percent.mt)
	}

	print( "load_data_PCA" )
	seurat_list = load_cellranger_result( sample_name, sample_dir_list, opt$o, max.percent.mt, min.nFeature_RNA, max.nFeature_RNA, min.features, file_name, opt$d )
	for (i in 1:length(seurat_list)){
		qc_metrics_plots( seurat_list[[i]], opt$o, sample_name[i] )
		seurat_list[[i]] <- filter_seurat_object( seurat_list[[i]], max.percent.mt, min.nFeature_RNA, max.nFeature_RNA  )
		variable_features_plots( seurat_list[[i]], opt$o, sample_name[i]  )
	}
	if ( length(seurat_list) == 1 ){
		seurat_object = seurat_list[[1]]
	}else{
		seurat_list.anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:10)
		seurat_object <- IntegrateData(anchorset = seurat_list.anchors, dims = 1:10)
	}
	seurat_object <- CellCycleScoring(object = seurat_object, s.features  = s.genes, g2m.features  = g2m.genes)
	seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
	all.genes <- rownames(seurat_object[["RNA"]])
	seurat_object <- ScaleData(seurat_object, features = all.genes, vars.to.regress = c("percent.mt","CC.Difference"))
	seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = TRUE)
	visualize_PCA( seurat_object, opt$o, file_name  )
	saveRDS( seurat_object, file = paste( opt$o, "/", file_name, "_load_data_PCA.rds" ,sep = ""))
	print( "请根据ElbowPlot图选择合适的PC参数并在ini文件中设置dim参数。可以参考示例文档" )
}

if ( Sample.ini$Analysis$analysis == "cluster" ){ #根据选择的PC数进行聚类，并绘制不同分辨率下的聚类图
	if ( is.null(opt$r) ){
		rds_file <- paste( opt$o, "/", file_name, "_load_data_PCA.rds" ,sep = "")
	}else{
		rds_file <- opt$r
	}
	if (file.exists(rds_file)){
		print( "cluster" )
	}else{
		print("请指定rds文件路径")
		q()
	}

	if ( is.null(Sample.ini$Cluster_Parameter$features) || is.null(Sample.ini$Cluster_Parameter$dim) || is.null(Sample.ini$Cluster_Parameter$resolution) ){ #读取配置文件中的参数
		print("请在配置文件中配置[Cluster_Parameter]:features、dim、resolution三个参数，参见示例文档")
		usage()
		q()
	}else{
		features = strsplit(Sample.ini$Cluster_Parameter$features, ",")[[1]]
		res = as.numeric(strsplit(Sample.ini$Cluster_Parameter$resolution, ",")[[1]])
		dims = as.numeric(strsplit(Sample.ini$Cluster_Parameter$dim, ",")[[1]])
	}
	file_name_temp <- file_name
	for (i in 1:length(dims)){
		seurat_object <- readRDS( file = rds_file )
		dim = dims[i]
		file_name = paste( file_name_temp, "_Dim", dim, sep = "" )
		seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:dim)
		seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:dim)
		seurat_object <- FindClusters(seurat_object, resolution = res)
		sample_cluster_DimPlot( seurat_object, opt$o, file_name, res, sample_name )
		DefaultAssay(seurat_object) <- "RNA"
		nk.markers = FindAllMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, grouping.var = "stim")
		write.table(nk.markers ,file = paste( opt$o, "/", file_name, "_all_marker.txt" ,sep = "" ), quote = FALSE,row.names = TRUE, sep = "\t")
		features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
		saveRDS( seurat_object, file = paste( opt$o, "/", file_name, "_cluster.rds" ,sep = ""))
		features_plot( seurat_object, nk.markers, features, opt$o, file_name )
	}
	file_name <- file_name_temp
}

if ( Sample.ini$Analysis$analysis == "mark_mutation_and_celltype" ){ #点亮突变信息及标记细胞类型
	if (is.null( Sample.ini$Mark_Mutation_And_Celltype_Parameter$dim ) || is.null( Sample.ini$Mark_Mutation_And_Celltype_Parameter$resolution ) || is.null( Sample.ini$CellType )){
		print("请在配置文件中配置[Mark_Mutation_And_Celltype_Parameter]:dim、resolution及[CellType]参数，参见示例文档")
		usage()
		q()
	}else{
		res = as.numeric(Sample.ini$Mark_Mutation_And_Celltype_Parameter$resolution)
		dim = as.numeric(Sample.ini$Mark_Mutation_And_Celltype_Parameter$dim)
	}
	file_name_temp <- file_name
	file_name = paste( file_name_temp, "_Dim", dim, sep = "" )
	if ( is.null(opt$r) ){
		rds_file <- paste( opt$o, "/", file_name, "_cluster.rds" ,sep = "")
	}else{
		rds_file <- opt$r
	}
	if (file.exists(rds_file)){
		seurat_object <- readRDS( file = rds_file )
	}else{
		print("请指定rds文件路径")
		q()
	}
	
	file_name = paste( file_name, "_Res", res, sep = "" )
	eval(parse(text = paste("Idents(object = seurat_object) <- seurat_object@meta.data$integrated_snn_res.",res ,sep = "")))
	seurat_object$seurat_clusters <- Idents(object = seurat_object)
	new.cluster.ids <- Cluster_cellType( Sample.ini$CellType, levels(seurat_object) )
	names(new.cluster.ids) <- levels(seurat_object)
	seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
	plots = DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.3, label.size = 1.5) + NoLegend()
	ggsave( paste( opt$o, "/", file_name, "_celltype_DimPlot.pdf" ,sep = "" ), plots, width = 10, limitsize = FALSE , height = 10)
	summary_matrix <- cluster_calculation( seurat_object, paste( opt$o, "/", file_name, "_summary_matrix.txt" ,sep = "" ), sample_name)
	summary_matrix[ is.na(summary_matrix) ] =0
	cell_no_barplot( summary_matrix, opt$o, file_name )
	percent_matrix <- summary_matrix/colSums(summary_matrix)
	cell_no_barplot( percent_matrix, opt$o, paste(file_name, "_percent",sep = "") )
	features = strsplit(Sample.ini$Cluster_Parameter$features, ",")[[1]]
	features <- features[features %in% rownames(seurat_object[["RNA"]]@data)]
	plots <- DotPlot(object = seurat_object, features = features) + RotatedAxis()
	ggsave( paste( opt$o, "/", file_name, "_DotPlot_marker.pdf" ,sep = "" ), plots, width=60,height=25,limitsize = FALSE )
	
	for ( sample_gene in  names(Sample.ini$Mutation)){
		dir <- Sample.ini$Mutation[[sample_gene]]
		gene <- strsplit( sample_gene, "_" )[[1]][2]
		gene_mutation = read.table( dir )
		colnames_result<-paste( "finalCellMutaionStatus", gene, sep="_")
		if ( colnames_result %in% names(seurat_object@meta.data)){
			result = seurat_object@meta.data[[colnames_result]]
		}else{
			result = as.data.frame(rep( "NoCoverage",length( seurat_object$stim ) ))
			rownames(result) = names(seurat_object$stim)
		}
		test <- as.matrix(subset(gene_mutation,select=c("finalCellMutaionStatus")))
		result <- as.matrix( result )
		row_names <- intersect(rownames(test),rownames(result))
		result[row_names,1] <- test[row_names,"finalCellMutaionStatus"]
		
		colnames(result) <- colnames_result
		seurat_object[[colnames_result]] <- result
	}
	
	mutation_marker_plot( Sample.ini$Mutation, seurat_object, paste( opt$o, "/", file_name, "_mutation_marker.png" ,sep = "" ), res )
	mutation_marker_plot_NoLegend( Sample.ini$Mutation, seurat_object, paste( opt$o, "/", file_name, "_mutation_marker_NoLegend.png" ,sep = "" ), res )
	plots = DimPlot(seurat_object, reduction = "umap", split.by = "stim",pt.size = 0.1)
	ggsave( paste( opt$o, "/", file_name, "_sample_DimPlot.pdf" ,sep = "" ), plots, width = 10, limitsize = FALSE, height = 10 )
	plots = DimPlot(seurat_object, reduction = "umap", split.by = "Phase",pt.size = 0.1)
	ggsave( paste( opt$o, "/", file_name, "_Phase_DimPlot.pdf" ,sep = "" ), plots, width = 33, limitsize = FALSE, height = 10 )
	saveRDS( seurat_object, file = paste( opt$o, "/", file_name, "_mark_mutation_and_celltype.rds", sep = ""))
}

#从数据中画gene pattern
library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
usage <- function(){ #使用说明
	cat("\n使用说明\n")
	print_help(opt_parser)
}

option_list = list(
	make_option( c("-i", "--in"), type="character", default=NULL, help="[Required] 配置文件，输入需要绘制的gene list", metavar = "character"),
	make_option( c("-f", "--foldChange"), type="character", default=2, help="[Required] 判断表达上调下降时的阈值,要求大于1", metavar="character"),
	make_option( c("-o", "--outdir"), type="character", default=NULL, help="[Required] 结果输出目录", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if ( is.null(opt$i) ){ 
	cat("\n配置文档格式：\n")
	cat("
[GeneList]
gene_list1=APAF1,TP53AIP1,BAX,PERP,FAS,PIDD,DR5,PIG3,PUMA,NOXA,SIVA
gene_list2=APAF1,TP53AIP1,BAX,PERP,FAS,PIDD,DR5,PIG3,PUMA,NOXA,SIVA
[Sample]
P101=/pnas/wangqf_group/suyx/PMO/test/genePattern/P101SR_26SampleAlignmentSeurat.rds
[Time]
D0=P101SR00-PB
D1=P101SR01-PB
D4=P101SR04-PB
D10=P101SR10-PB
")
	print_help(opt_parser)
	q()
}
if ( is.null(opt$o)){
	cat( "\n请制定输出目录\n" )
	print_help(opt_parser)
	q()
}

library( dplyr )
library( Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library( ggplot2 )
library( cowplot )
opt$f = as.numeric(opt$f)
Sample.ini <- read.ini(opt$i)
sample_name = names(Sample.ini$Sample)
gene_list_name = names(Sample.ini$GeneList)
for (i in 1:length(sample_name) ){
	sample = sample_name[i]
	rdsfile = Sample.ini$Sample[i][[1]]
	seurat_object = readRDS( file = rdsfile )
	express_matrix = as.matrix(seurat_object[["RNA"]]@data)
	clusters = unique(as.character(seurat_object@active.ident))
	cells = names( seurat_object@active.ident )
	cells_cluster = as.character(seurat_object@active.ident)
	samples = seurat_object@meta.data$stim
	names( cells_cluster ) = cells
	names(samples) = cells

	for (j in 1:length( gene_list_name )){
		gene_list_row = strsplit( Sample.ini$GeneList[i][[1]], "," )[[1]]
		gene_list = gene_list_row[ which( gene_list_row %in% rownames( express_matrix ) ) ]
		express_matrix_plot = express_matrix[ gene_list, ]
		plot_col_num = length( clusters )
		plot_row_num = 4
		pdf( paste(opt$o, "/", sample, "_", gene_list_name, ".pdf", sep = ""), width = plot_row_num*5, height = plot_col_num*5 )
		opar<-par(no.readonly=TRUE)
		par(mfrow=c(plot_col_num,plot_row_num))
		for ( cluster in clusters ){
			cluster_express_matrix = express_matrix_plot[,names( which(cells_cluster == cluster) )]
			samples_cluster = samples[ names( which(cells_cluster == cluster) ) ]
			Time = names( Sample.ini$Time )
			up = down = undulation = flat = matrix( rep(0, length(Time)), 1, length(Time))
			for (n in 1:dim(cluster_express_matrix)[1]){
				temp = c()
				for (m in 1:length(Time)){
					sample_time_express = cluster_express_matrix[ n, names(samples_cluster[which(samples_cluster %in% strsplit(Sample.ini$Time[m][[1]], ","))])]
					sample_time_express = sample_time_express[ which(sample_time_express > 0) ]
					if ( sum(sample_time_express) == 0 ){
						temp = c(temp,0.0001)
					}else{
						temp = c(temp,mean(sample_time_express))
					}
					
				}
				status = "flat"
				for ( x in 2:length( temp ) ){
					if ( temp[x]/temp[x-1] >= opt$f ){
						if ( status == "flat" || status == "up"){
							status = "up"
						}else{
							status = "undulation"
						}
					}
					if ( temp[x]/temp[x-1] <= 1/opt$f ){
						if ( status == "flat" || status == "down"){
							status = "down"
						}else{
							status = "undulation"
						}
					}
				}
				if ( status == "up" ){ up = rbind(up,temp) }
				if ( status == "down" ){ down = rbind(down,temp) }
				if ( status == "undulation" ){ undulation = rbind(undulation,temp) }
				if ( status == "flat" ){ flat = rbind(flat,temp) }
			}
			plot( 1:dim(up)[2], up[1,], main = cluster, ylab = "up", type = "l", ylim = c(0,1),lwd=0, yaxs="i",xaxs="i",xaxt="n")
			axis(1,1:dim(flat)[2],labels=Time)
			if (dim(up)[1] > 1){
				for ( y in 2:dim(up)[1]){
					lines( 1:dim(up)[2], up[y,] )
				}
			}
			
			plot( 1:dim(down)[2], down[1,], main = cluster, ylab = "down", type = "l", ylim = c(0,1),lwd=0, yaxs="i",xaxs="i",xaxt="n")
			axis(1,1:dim(flat)[2],labels=Time)
			if (dim(down)[1] > 1){
				for ( y in 2:dim(down)[1]){
					lines( 1:dim(down)[2], down[y,] )
				}
			}
			
			plot( 1:dim(undulation)[2], undulation[1,], main = cluster, ylab = "undulation", type = "l", ylim = c(0,1),lwd=0, yaxs="i",xaxs="i",xaxt="n")
			axis(1,1:dim(flat)[2],labels=Time)
			if (dim(undulation)[1] > 1){
				for ( y in 2:dim(undulation)[1]){
					lines( 1:dim(undulation)[2], undulation[y,] )
				}
			}
			
			plot( 1:dim(flat)[2], flat[1,], main = cluster, ylab = "flat", type = "l", ylim = c(0,1) ,lwd=0, yaxs="i",xaxs="i",xaxt="n")
			axis(1,1:dim(flat)[2],labels=Time)
			if (dim(flat)[1] > 1){
				for ( y in 2:dim(flat)[1]){
					lines( 1:dim(flat)[2], flat[y,] )
				}
			}
		}
		par(opar)
		dev.off()
	}
}


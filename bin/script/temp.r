##单独运行
library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)
dims=c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
res=c(0.5,0.8,1.0,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8)
features<-c("CEACAM1","CEACAM6","CEACAM3","MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","CD38","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CCR7","CD4","FOXP3","SELL","S100A4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","LILRA4","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")

width = 6
height = 6
sample_temp = "P106"
sample = sample_temp
outdir = paste( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/", sample, sep = "")
sample_name = 1:2
for ( dim in dims ){
	sample = paste( sample_temp, "_Dim", dim, sep = "" )
	seurat_object_file = paste( outdir, "/", sample, "_cluster.rds", sep = "")
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


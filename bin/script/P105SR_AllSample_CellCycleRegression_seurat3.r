#==========================================================================
#
# Goal：D0 ~ D26 alignment(P105SR,PB and BM ,6 samples，P105SR_AllSample_CellCyleRegression)
# By: 	Jiangst
# Date: 2019-04-29
# Software: Seurat
#==========================================================================

rm(list = ls())
library(Seurat,lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
packageVersion("Seurat")  #'3.0.0.9000'
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(cowplot)

wd<-"/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P105SR_AllSample_CellCyleRegression"
#dim<-args[2]
setwd(wd)
patientID<-"P105SR"
sampleName_temp<-dir(path="/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result",pattern="P105SR*")
sampleName<-substr(sampleName_temp,1,11)
sampleName

infor2count<-c("rawCellNo","remainCellNo","filtered","percent.mito>=0.15","nGene<=200")
count<-matrix(0,nrow=length(sampleName),ncol=length(infor2count))
rownames(count)<-sampleName
colnames(count)<-infor2count

sampleList_seurat=c()
for (i in 1:length(sampleName)) {
	a<-paste("Sample",i,sep="")
	sampleList_seurat<-c(sampleList_seurat,a)
}
sampleList_seurat

for  (j in 2:length(sampleName)) {
	Path<-paste("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/Link2AllSampleResult/ExpressionMatrix/",sampleName[j],sep='')
	list.files(Path)
	Sample_10x<-Read10X(data.dir = Path)
	Sample<- CreateSeuratObject(counts = Sample_10x,min.cells=3,min.features=200,project = sampleName[j])
    Sample<- RenameCells(Sample,add.cell.id=sampleName[j])  #check CB:colnames(x = Sample)[1:3]
	Sample$stim <- sampleName[j]
    Sample[["percent.mt"]] <- PercentageFeatureSet(object = Sample, pattern = "^MT-")
	##统计过滤前后的信息
	count[j,1]<-dim(Sample)[2]
	count[j,2]<-dim(subset(Sample@meta.data,nFeature_RNA>200 & percent.mt<15))[1]
	count[j,3]<-dim(subset(Sample@meta.data,nFeature_RNA<=200))[1]+dim(subset(Sample@meta.data,percent.mt>=15))[1]-dim(subset(Sample@meta.data,nFeature_RNA<=200 & percent.mt>=15))[1]
	count[j,4]<-dim(subset(Sample@meta.data,percent.mt>=15))[1]
	count[j,5]<-dim(subset(Sample@meta.data,nFeature_RNA<=200))[1]
    Sample <- subset(x = Sample, subset = nFeature_RNA > 200  & percent.mt < 15)   # nFeature_RNA < 2500
	#Seurat Normalizing the data
    Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize", scale.factor = 10000,verbose=TRUE)
    Sample <- FindVariableFeatures(object = Sample, selection.method = "vst", nfeatures = 2000)  # How to choose top variable features:vst, mean.var.plot,dispersion
	assign(sampleList_seurat[j],Sample)
	rm(Sample_10x)
}

count
write.table(count,paste(patientID,"_SampleFilteringCount.txt",sep=""),sep="\t",quote = FALSE)


#Perform integration
SampleAlignmentSeurat.list <- list(get(sampleList_seurat[1]),get(sampleList_seurat[2]),get(sampleList_seurat[3]),get(sampleList_seurat[4]),get(sampleList_seurat[5]),get(sampleList_seurat[6]))
SampleAlignmentSeurat.anchors<- FindIntegrationAnchors(object.list = SampleAlignmentSeurat.list, dims = 1:25)
SampleAlignmentSeurat <- IntegrateData(anchorset = SampleAlignmentSeurat.anchors, dims = 1:25)

#Perform an integrated analysis
DefaultAssay(object = SampleAlignmentSeurat) <- "integrated"
all.genes <- rownames(x = SampleAlignmentSeurat) #对所有gene标准化，doheatmap会需要

#calculate the cell cycle score
cc.genes <- readLines(con = "/asnas/wangqf_group/jiangsht/Test4Seurat/Src/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
SampleAlignmentSeurat <- CellCycleScoring(object = SampleAlignmentSeurat, s.features  = s.genes, g2m.features  = g2m.genes)
head(SampleAlignmentSeurat[[]])
SampleAlignmentSeurat$CC.Difference <- SampleAlignmentSeurat$S.Score - SampleAlignmentSeurat$G2M.Score
SampleAlignmentSeurat <- ScaleData(object = SampleAlignmentSeurat, features = all.genes, vars.to.regress = c("percent.mt","CC.Difference"))
SampleAlignmentSeurat <- RunPCA(object = SampleAlignmentSeurat, npcs = 30, verbose = TRUE)

#选择PC
pdf(file = paste(patientID,"_SampleAlignmentSeurat_PCDetermination_heatmap.pdf",sep=""),width=14,height = 30)
DimHeatmap(object = SampleAlignmentSeurat,  dims = 1:30, cells = 1000, balanced = TRUE)
dev.off()

pdf(file = paste(patientID,"_SampleAlignmentSeurat_PCDetermination_DimLoadings.pdf",sep=""),width=14,height = 50)
VizDimLoadings(object = SampleAlignmentSeurat, dims = 1:30)
dev.off()

pdf(file = paste(patientID,"_SampleAlignmentSeurat_CCADetermination_Elbow.pdf",sep=""),width=14,height = 8)
ElbowPlot(SampleAlignmentSeurat, ndims = 30)
dev.off()

saveRDS(SampleAlignmentSeurat, file = paste(patientID,"_SampleAlignment.rds",sep=""))



####=================
#19-23
dim=20
res=c(1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.6,3.8)
wd<-"/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample_seurat3.0/P105SR_AllSample_CellCyleRegression"
setwd(wd)
patientID<-"P105SR"
sampleName_temp<-dir(path="/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result",pattern="P105SR*")
sampleName<-substr(sampleName_temp,1,11)
sampleName

rdsfile<-paste(wd,"/",patientID,"_SampleAlignment.rds",sep="")
SampleAlignmentSeurat <- readRDS(file = rdsfile)
SampleAlignmentSeurat <- RunUMAP(object = SampleAlignmentSeurat,reduction = "pca", dims = 1:dim)
SampleAlignmentSeurat <- FindNeighbors(object = SampleAlignmentSeurat, reduction = "pca", dims = 1:dim)
SampleAlignmentSeurat <- FindClusters(SampleAlignmentSeurat, resolution = res)

#=================================mark the mutation=================================
#add PTPN11
gene="PTPN11"
mutationInfo<-read.table(file=paste("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample/Test/Test4CountMutation/ExtractCBFromScRNA/PysamExtractReads/",patientID,"MutationDetection_200gene_15mit_withoutTrim_",gene,".txt",sep=""),header=TRUE,sep="\t")
test<-subset(mutationInfo,select=c("finalCellMutaionStatus"))
colnames(test)<-paste(colnames(test),gene,sep="_")
SampleAlignmentSeurat$finalCellMutaionStatus_PTPN11<- test

#add FLT3
gene="FLT3"
mutationInfo<-read.table(file=paste("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample/Test/Test4CountMutation/ExtractCBFromScRNA/PysamExtractReads/",patientID,"MutationDetection_200gene_15mit_withoutTrim_",gene,".txt",sep=""),header=TRUE,sep="\t")
test<-subset(mutationInfo,select=c("finalCellMutaionStatus"))
colnames(test)<-paste(colnames(test),gene,sep="_")
SampleAlignmentSeurat$finalCellMutaionStatus_FLT3<- test

#=================================Visualization=================================
#show the different resolution
nrow= ceiling(length(res)/4)
p<-c()
for (i in 1:length(res)) {
	a<-paste("p",i,sep="")
	p<-c(p,a)
}
for (i  in 10:(dim(SampleAlignmentSeurat@meta.data)[2]-5)) {	
	temp <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 0.5, label = TRUE, label.size=6,group.by = colnames(SampleAlignmentSeurat@meta.data[i]),order=c(levels(factor(SampleAlignmentSeurat@meta.data[i],levels=c((dim(unique(SampleAlignmentSeurat@meta.data[i]))[1]-1):0))))) +NoLegend() + labs(title=colnames(SampleAlignmentSeurat@meta.data[i]))
	assign(p[i-9],temp)
}

png(file = paste(patientID,"_",dim,"DimPlot_differentRes.png",sep=""),width=2400,height = 800*nrow)
plot_grid(get(p[1]),get(p[2]),get(p[3]),get(p[4]),get(p[5]),get(p[6]),get(p[7]),get(p[8]), get(p[9]),get(p[10]),get(p[11]),get(p[12]),get(p[13]),get(p[14]),get(p[15]),get(p[16]),get(p[17]),get(p[18]),get(p[19]),get(p[20]),nrow=nrow,ncol=4)
dev.off()

#Visualization 
res=2
Idents(object = SampleAlignmentSeurat) <- SampleAlignmentSeurat@meta.data$integrated_snn_res.2 #change the ident
SampleAlignmentSeurat$seurat_clusters <- Idents(object = SampleAlignmentSeurat)

png(file = paste(patientID,"_",dim,"UMAP_res",res,".png",sep=""),width=1400,height = 2400)
#pdf(file = paste(dim,"Samplemerge_cellCycle_TSNEPlot_res",res,".pdf",sep=""),width=25,height = 8)
p1 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "stim", order=names(sort(table(SampleAlignmentSeurat@meta.data$stim))), pt.size = 1) #different sample
p2 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 1, label = TRUE, label.size=6) # +NoLegend() #different cluster
p3 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "Phase",pt.size = 1) #, order=c("G1","S","G2M")
p4 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "finalCellMutaionStatus_PTPN11", order=c("Mut","WT","NoCoverage"), cols=c("#CCCCCC","#FFFF00","#CC0033"),pt.size = 2.5)
p5 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "finalCellMutaionStatus_FLT3", order=c("Mut","WT","NoCoverage"), cols=c("#CCCCCC","#FFFF00","#CC0033"),pt.size = 2.5)
plot_grid(p1, p2,p3,p4,p5,ncol=2)
dev.off()

#split by sample
png(file = paste(patientID,"_",dim,"UMAP_splitBySample_res",res,".png",sep=""),width=2200,height = 2200)
DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 2.5,split.by="stim", label = TRUE, label.size=6)+NoLegend() # +NoLegend() #different cluster
dev.off()

#split by cell cycle
png(file = paste(patientID,"_",dim,"UMAP_splitByCellCycle_res",res,".png",sep=""),width=2100,height = 800)
DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 1.5,split.by="Phase", label = TRUE)+NoLegend() # +NoLegend() #different cluster
dev.off()

#marker plot
DefaultAssay(object = SampleAlignmentSeurat) <- "RNA"
features<-c("CEACAM1","CEACAM6","CEACAM3","MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","CD38","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CCR7","CD4","FOXP3","SELL","S100A4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","LILRA4","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")
features[!features %in% rownames(SampleAlignmentSeurat[["RNA"]]@data)]
features<-features[features %in% rownames(SampleAlignmentSeurat[["RNA"]]@data)]
nrow= ceiling(length(features)/4)

png(file = paste(patientID,"_",dim,"Samplemerge_cluster_marker_res",res,".png",sep=""),width=2400,height = 800*nrow,pointsize=15)
FeaturePlot(object = SampleAlignmentSeurat, features = features,pt.size = 1,min.cutoff="q9") #, cols = c("#CCCCCC", "#FF0000")
dev.off()	
png(file = paste(patientID,"_",dim,"Samplemerge_VlnPlot_marker_res",res,".png",sep=""),width=3360,height =320*nrow ,pointsize=10)
VlnPlot(object = SampleAlignmentSeurat, pt.size = 0.5,features =features,ncol=4)
dev.off()
p <- DotPlot(object = SampleAlignmentSeurat, features = features) + RotatedAxis()
ggsave(p,filename=paste(patientID,"_",dim,"Samplemerge_DotPlot_marker_res",res,".pdf",sep=""),width=60,height=25,limitsize = FALSE)

saveRDS(SampleAlignmentSeurat, file = paste(patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))

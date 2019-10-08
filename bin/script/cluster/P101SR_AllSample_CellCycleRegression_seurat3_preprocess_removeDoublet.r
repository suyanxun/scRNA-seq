#==========================================================================
#
# Goal：D0 ~ D26 alignment,5 samples(P101SR) for diff CC and res; remove doublet
# By: 	Jiangst
# Date: 2019-06-12
# Update:Seurat 3.0; preprocess bam
#==========================================================================

#加载需要的library
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

wd <-"/asnas/wangqf_group/zhouzihan/scRNA/Seurat/CellCycleRegressionRomoveDB"
setwd(wd)
patientID <- "P108SR"
sampleName <- dir(path="/asnas/wangqf_group/zhouzihan/scRNA/Seurat/ExpressionMatrix2",pattern="P108SR*")

infor2count <- c("rawCellNo","remainCellNo","totalFiltered","percent.mito>=0.15","nGene<=200","predictDoubletNo","predictDoubletRatio")
count <- matrix(0,nrow=length(sampleName),ncol=length(infor2count))
rownames(count) <- sampleName
colnames(count) <- infor2count

SampleAlignmentSeurat.list =c()
allDoublet <- c()
for  (j in 1:length(sampleName)) {
	Path<-paste("/asnas/wangqf_group/zhouzihan/scRNA/Seurat/ExpressionMatrix2/",sampleName[j],"/",sep='')
	list.files(Path)
	Sample_10x<-Read10X(data.dir = Path) #Read10X将3个矩阵读成稀疏矩阵，CellRanger3.0后输入格式必须为.gz
	Sample<- CreateSeuratObject(counts = Sample_10x,min.cells=3,min.features=200,project = sampleName[j]) #创建SeuratObject，此时初步过滤细胞
    Sample<- RenameCells(Sample,add.cell.id=sampleName[j])  #check CB:colnames(x = Sample)[1:3],在barcode前加上样本名称
    count[j,1]<-dim(Sample)[2]  #将现在样本的细胞数量记录到count矩阵中
    #extract doublet
    doubletStatus <- read.table(paste("/asnas/wangqf_group/zhouzihan/scRNA/scrublet/Result/",sampleName[j],"/",sampleName[j],".predictDoublet_scrublet.txt",sep = ""))
    cellBarcode <- read.table(paste("/asnas/wangqf_group/zhouzihan/scRNA/Seurat/ExpressionMatrix2/",sampleName[j],"/barcodes.tsv", sep = ""))
    cellBarcode <- paste(sampleName[j], substr(cellBarcode[,1],1,16), sep = "_") #在barcode前加上样本名称
    doublet <- data.frame(cellBarcode = cellBarcode, doubletStatus_Scrublt = doubletStatus)
    colnames(doublet) <- c("cellBarcode","doubletStatus_Scrublt")
    rownames(doublet) <- cellBarcode
    temp <- doublet[which(doublet$doubletStatus_Scrublt == "True"),"cellBarcode"] #将预测为doublet的细胞存储进temp中
    temp <- as.character(temp)
    allDoublet <- c(allDoublet,temp)
    removeDoublet <- setdiff(colnames(Sample), temp)
    count[j,6] <- length(colnames(Sample)) - length(removeDoublet)
    #remove doublet
    Sample <- subset(Sample, cells = removeDoublet)
	Sample$stim <- sampleName[j]
    Sample[["percent.mt"]] <- PercentageFeatureSet(object = Sample, pattern = "^MT-") #在meta.data中加入percent.mt，记录MT基因占总基因的比例
    plots = VlnPlot(object = Sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave( paste( sampleName[j],"_All_features.pdf" ,sep = "" ), plots,width=15,height=8)
	plot1 <- FeatureScatter(object = Sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
	plot2 <- FeatureScatter(object = Sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	plots = CombinePlots(plots = list(plot1, plot2))
	ggsave( paste( sampleName[j], "_CombinePlots.pdf" ,sep = "" ), plots ,width=15,height=8)
	count[j,2] <- dim(subset(Sample@meta.data,nFeature_RNA>200 & percent.mt<15))[1] 
	count[j,3] <- dim(subset(Sample@meta.data,nFeature_RNA<=200))[1]+dim(subset(Sample@meta.data,percent.mt>=15))[1]-dim(subset(Sample@meta.data,nFeature_RNA<=200 & percent.mt>=15))[1] + count[j,6]
	count[j,4] <- dim(subset(Sample@meta.data,percent.mt>=15))[1]
	count[j,5] <- dim(subset(Sample@meta.data,nFeature_RNA<=200))[1]   
    Sample <- subset(x = Sample, subset = nFeature_RNA > 200  & percent.mt < 15)   # 过滤掉低质量细胞
	#Seurat Normalizing the data
    Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize", scale.factor = 10000,verbose=TRUE)
    Sample <- FindVariableFeatures(object = Sample, selection.method = "vst", nfeatures = 2000)  # How to choose top variable features:vst, mean.var.plot,dispersion
	SampleAlignmentSeurat.list =c(SampleAlignmentSeurat.list,Sample)
	rm(Sample_10x)
}

count
count[,7] <- count[,6] / count[,1] *100
write.table(count,paste(patientID,"_SampleFilteringCount.txt",sep=""),sep="\t",quote = FALSE)
write.table(allDoublet, file = paste(patientID,"_doubletScrublt.txt",sep=""),sep="\t",quote = FALSE, row.names = FALSE , col.names = FALSE)

#================================Perform integration==============
SampleAlignmentSeurat.anchors<- FindIntegrationAnchors(object.list = SampleAlignmentSeurat.list, dims = 1:30)
SampleAlignmentSeurat <- IntegrateData(anchorset = SampleAlignmentSeurat.anchors, dims = 1:30)

#Perform an integrated analysis
all.genes <- rownames(SampleAlignmentSeurat[["RNA"]]) #对所有gene标准化，doheatmap会需要
DefaultAssay(object = SampleAlignmentSeurat) <- "integrated"

#calculate the cell cycle score
cc.genes <- readLines(con = "/asnas/wangqf_group/jiangsht/Test4Seurat/Src/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
SampleAlignmentSeurat <- CellCycleScoring(object = SampleAlignmentSeurat, s.features  = s.genes, g2m.features  = g2m.genes)
head(SampleAlignmentSeurat[[]])
SampleAlignmentSeurat$CC.Difference <- SampleAlignmentSeurat$S.Score - SampleAlignmentSeurat$G2M.Score
SampleAlignmentSeurat <- ScaleData(object = SampleAlignmentSeurat, features = all.genes, vars.to.regress = c("percent.mt","CC.Difference"))
SampleAlignmentSeurat <- RunPCA(object = SampleAlignmentSeurat, npcs = 50, verbose = TRUE)

#选择PC
pdf(file = paste(patientID,"_SampleAlignmentSeurat_PCDetermination_heatmap.pdf",sep=""),width=14,height = 30)
DimHeatmap(object = SampleAlignmentSeurat,  dims = 1:50, cells = 1000, balanced = TRUE)
dev.off()

pdf(file = paste(patientID,"_SampleAlignmentSeurat_PCDetermination_DimLoadings.pdf",sep=""),width=14,height = 50)
VizDimLoadings(object = SampleAlignmentSeurat, dims = 1:50)
dev.off()

pdf(file = paste(patientID,"_SampleAlignmentSeurat_PCDetermination_Elbow.pdf",sep=""),width=14,height = 8)
ElbowPlot(SampleAlignmentSeurat, ndims = 50)
dev.off()

pdf(file = paste(patientID,"_PCDetermination_JackStraw.pdf",sep=""),width=14,height = 8)
SampleAlignmentSeurat <- JackStraw(SampleAlignmentSeurat, num.replicate = 100,dims=50)
SampleAlignmentSeurat <- ScoreJackStraw(SampleAlignmentSeurat, dims = 1:50)
JackStrawPlot(SampleAlignmentSeurat, dims = 1:50)
dev.off()

saveRDS(SampleAlignmentSeurat, file = paste(patientID,"_SampleAlignment.rds",sep=""))

#========================================================
dim=21
res=2.4
wd <-"/asnas/wangqf_group/zhouzihan/scRNA/Seurat/NewIdent"
setwd(wd)
patientID <- "P108SR"
sampleName <- dir(path="/asnas/wangqf_group/zhouzihan/scRNA/Seurat/ExpressionMatrix2",pattern="P108SR*")
SampleAlignmentSeurat <- readRDS(file = paste('/asnas/wangqf_group/zhouzihan/scRNA/Seurat/DiffPC/', patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))
Idents(object=SampleAlignmentSeurat) <- paste("integrated_snn_res.",res,sep="") #??

#==========================change ident===================================================
new.ident <- c("stem_cell","prog","prog","stem_cell","prog","unknown1","NK","Tm","ery","Tm","CD16-Mono","CD4T","prog","stem_cell","stem_cell","prog","CD8T","CD16-Mono","CD8T","NK","stem_cell","Tm","stem_cell","prog","CD8T","CD16-Mono","unknown2","gran","ery","ery","CD16-Mono","CD16+Mono","unknown1","ery","ery","ery","NK","plasmaB","mega","CD1C+DC","normalB","unknown2","prog","pDC","naiveT","naiveT","unknown3")
new.ident <- c("prog","Tm","CD4T","Eryth","CD14mono_1","prog","granu","CD4T","stemCell_2","unknown","CD14mono_2","stemCell_1","GMP","Tm","prog","Tm","NK","Tm","CD14mono_2","plasmaB","MK","CD14mono_1","CD8T","Tm","CD16mono","stemCell_1","CD8T","GMP","granu","granu","CD8T","GMP","plasmaB","stemCell_3","pDC","normalB")
length(new.ident)  ##30
for (i in 0:(length(new.ident)-1)) {
    subCell<-WhichCells(SampleAlignmentSeurat,idents=i)
    SampleAlignmentSeurat <- SetIdent(object = SampleAlignmentSeurat, cells = subCell, value = new.ident[i+1])
}
SampleAlignmentSeurat$orig.ident<-Idents(object = SampleAlignmentSeurat)

png(file = paste(patientID,"_",dim,"UMAP_ChangeIdent_withoutIdent_beforeMerge_res",res,".png",sep=""),width=2000,height = 2400)
#pdf(file = paste(dim,"Samplemerge_cellCycle_TSNEPlot_res",res,".pdf",sep=""),width=25,height = 8)
p1 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "stim", order=names(sort(table(SampleAlignmentSeurat@meta.data$stim))), pt.size = 1) +NoLegend()  #different sample
p2 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 1, label = FALSE, label.size=8) +NoLegend()
p3 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "Phase",pt.size = 1) +NoLegend() 
#p4 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = paste("finalCellMutaionStatus_",gene[1],sep=""), order=c(paste(gene[1],c("Mut","WT","further2confirm"),sep="_"),"NoCoverage"), cols=c("#CCCCCC","#CCCCCC","#0000FF","#CC0033"),pt.size = 2.5)+NoLegend() 
#p5 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = paste("finalCellMutaionStatus_",gene[2],sep=""), order=c(paste(gene[2],c("Mut","WT","further2confirm"),sep="_"),"NoCoverage"), cols=c("#CCCCCC","#CCCCCC","#0000FF","#CC0033"),pt.size = 2.5)+NoLegend() 
plot_grid(p1, p2,p3,ncol=2) +NoLegend() 
dev.off()

png(file = paste(patientID,"_",dim,"UMAP_ChangeIdent_beforeMerge_res",res,".png",sep=""),width=2000,height = 2400)
p1 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "stim", order=names(sort(table(SampleAlignmentSeurat@meta.data$stim))), pt.size = 1) 
p2 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 1, label = TRUE, label.size=8) +NoLegend()  
p3 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = "Phase",pt.size = 1)   #, order=c("G1","S","G2M")
#p4 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = paste("finalCellMutaionStatus_",gene[1],sep=""),  order=c(paste(gene[1],c("Mut","WT","further2confirm"),sep="_"),"NoCoverage"), cols=c("#CCCCCC","#CCCCCC","#0000FF","#CC0033"),pt.size = 2.5) 
#p5 <- DimPlot(object = SampleAlignmentSeurat, reduction="umap",group.by = paste("finalCellMutaionStatus_",gene[2],sep=""),  order=c(paste(gene[2],c("Mut","WT","further2confirm"),sep="_"),"NoCoverage"), cols=c("#CCCCCC","#CCCCCC","#0000FF","#CC0033"),pt.size = 2.5) 
plot_grid(p1, p2,p3,ncol=2) 
dev.off()

#split by sample
png(file = paste(patientID,"_",dim,"UMAP_splitBySample_res",res,".png",sep=""),width=6600,height = 2200)
DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 2.5,split.by="stim", label = TRUE, label.size=6)+NoLegend() # +NoLegend() #different cluster
dev.off()

#split by sample and phase
png(file = paste(patientID,"_",dim,"UMAP_splitBySample_Phase_res",res,".png",sep=""),width=6600,height = 2200)
DimPlot(object = SampleAlignmentSeurat, reduction="umap",pt.size = 2.5,split.by="stim",group.by = "Phase", label = TRUE, label.size=6) # +NoLegend() #different cluster
dev.off()

DefaultAssay(object = SampleAlignmentSeurat) <- "RNA"
features<-c("CLEC12A","MECOM","CEACAM1","CEACAM6","CEACAM3","MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","CD38","THY1","ITGA6","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CCR7","CD4","FOXP3","SELL","S100A4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","LILRA4","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")
#features[!features %in% rownames(SampleAlignmentSeurat[["RNA"]]@data)]
features<-features[features %in% rownames(SampleAlignmentSeurat[["RNA"]]@data)]
nrow= ceiling(length(features)/4)

png(file = paste(patientID,"_",dim,"Samplemerge_cluster_marker_res",res,".png",sep=""),width=2400,height = 800*nrow,pointsize=15)
FeaturePlot(object = SampleAlignmentSeurat, features = features,pt.size = 1,min.cutoff="q9",order=TRUE) #, cols = c("#CCCCCC", "#FF0000")
dev.off()

png(file = paste(patientID,"_",dim,"Samplemerge_VlnPlot_marker_res",res,".png",sep=""),width=3360,height =320*nrow ,pointsize=10)
VlnPlot(object = SampleAlignmentSeurat, pt.size = 0.5,features =features,ncol=4, group.by = paste("integrated_snn_res.",res,sep=""))
dev.off()
p <- DotPlot(object = SampleAlignmentSeurat, features = features) + RotatedAxis()
ggsave(p,filename=paste(patientID,"_",dim,"Samplemerge_DotPlot_marker_res",res,".pdf",sep=""),width=60,height=25,limitsize = FALSE)

saveRDS(SampleAlignmentSeurat, file = paste(patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))

#=================================find the cluster marker in every cluster============================
DefaultAssay(object = SampleAlignmentSeurat) <- "integrated"
Idents(object=SampleAlignmentSeurat) <- paste("integrated_snn_res.",res,sep="") 
ClusterMarker <- FindAllMarkers(SampleAlignmentSeurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top50 <- ClusterMarker %>% group_by(cluster) %>% filter(avg_logFC > 0.5) %>% filter(p_val_adj < 0.01) %>% top_n(n = 50, wt = avg_logFC) 
top20 <- ClusterMarker %>% group_by(cluster) %>% filter(avg_logFC > 0.5) %>% filter(p_val_adj < 0.01) %>% top_n(n = 20, wt = avg_logFC) 
top10 <- ClusterMarker %>% group_by(cluster) %>% filter(avg_logFC > 0.5) %>% filter(p_val_adj < 0.01) %>% top_n(n = 10, wt = avg_logFC) 
write.table(top50,paste(patientID,"_",dim,"_Add_ClusterMarkerTop50_beforeMerge_",res,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
write.table(top20,paste(patientID,"_",dim,"_Add_ClusterMarkerTop20_beforeMerge_",res,".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

p <- DoHeatmap(object = SampleAlignmentSeurat, features  = top50$gene) 
q <- DoHeatmap(object = SampleAlignmentSeurat, features  = top20$gene) 
#ggsave( paste( patientID, "Samplemerge_heatmap_beforeMerge_",res,".png" ,sep = "" ), p ) #,width=15,height=8
png(file = paste(patientID,"_",dim,"_Add_heatmap_1_",res,".png",sep=""),width=1350,height = 2400)
plot(p)
dev.off()

png(file = paste(patientID,"_",dim,"_Add_heatmap_2_",res,".png",sep=""),width=1350,height = 2400)
plot(q)
dev.off()

saveRDS(SampleAlignmentSeurat, file = paste(patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))

#===============================================
#计算突变在的群体
SampleAlignmentSeurat$orig.ident<-Idents(object = SampleAlignmentSeurat)
#SampleAlignmentSeurat<-readRDS(file = paste(patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))

mutationNo <- table(SampleAlignmentSeurat@meta.data$orig.ident,SampleAlignmentSeurat@meta.data$finalCellMutaionStatus_CEBPA_ins)

#计算每个群体在每个样本中的突变数,用于确定肿瘤群体(D0>2)
mutationNoInClusterStim <-table(SampleAlignmentSeurat$orig.ident,SampleAlignmentSeurat$stim,SampleAlignmentSeurat$finalCellMutaionStatus_CEBPA_ins)
count <- mutationNoInClusterStim[,,1]

write.table(count, file=paste(patientID, 'MutationCount_', 'CEBPA_ins_', res, '.txt', sep=''), quote=F, sep='\t')

for(i in 1:length(gene)){
	mutationNo <- table(SampleAlignmentSeurat@meta.data$orig.ident,SampleAlignmentSeurat@meta.data$finalCellMutaionStatus_EP300)
	mutationNoInClusterStim <-table(SampleAlignmentSeurat$orig.ident,SampleAlignmentSeurat$stim,SampleAlignmentSeurat$finalCellMutaionStatus_EP300)
	count <- mutationNoInClusterStim[,,1]
	write.table(count, file=paste(patientID, 'MutationCount_', 'EP300_', res, '.txt', sep=''), quote=F, sep='\t')

}



#计算比例
clusterNo<-table(SampleAlignmentSeurat$orig.ident,SampleAlignmentSeurat$stim)
clusterRatio<-prop.table(clusterNo,2)*100

#计算比例（只看肿瘤群体）
#blast<-c("stemCell_1","pro_1","stemCell_2","cDC","pDC","stemCell_3","pro_2","pro_3")  #FLT3
blast<-c("stemCell_1","stemCell_2","stemCell_3","prog","GMP") #,,"granu"
colSums(clusterRatio[blast,])

#=================================ratio changes of each cluster=========================
#show the ratio changes of each cluster 
Freq <- table(SampleAlignmentSeurat$stim,SampleAlignmentSeurat$orig.ident,SampleAlignmentSeurat$Phase)
Freq<-as.data.frame(Freq)
colnames(Freq)<-c("SampleID","Cluster","Phase","Count")

Cluster_Percet = ddply(Freq,"SampleID",transform,percent = Count/sum(Count) * 100)
Cluster_Percet$SampleID<-factor(Cluster_Percet$SampleID,levels=sampleName)
Cluster_Percet$Phase<-factor(Cluster_Percet$Phase,levels=c("G2M","S","G1"))

colourCount = length(unique(Cluster_Percet$Cluster))
getPalette = colorRampPalette(colors = brewer.pal(9, "Set1"))

pdf(file = paste(patientID,"_",dim,"_Add_ClusterPercetAndCellCycleChanges_res",res,".pdf",sep=""),width=45,height = 25)
p1<-ggplot(data=Cluster_Percet, aes(x=SampleID, y=percent, fill=Cluster))+
  geom_bar(stat="identity")+
  #	geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
  labs(title=paste(dim,"clusterPercentChanges_res_",res,sep=""))+
  theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40),legend.position = "top")+
#  scale_fill_brewer(palette = "Set3",direction = -1)+
  scale_fill_manual(values = getPalette(colourCount)) +
  facet_grid(.~Cluster)
p2<-ggplot(data=Cluster_Percet, aes(x=SampleID, y=percent, fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  #	geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
  labs(title=paste(dim,"clusterCellCycleChanges_res_",res,sep=""))+
  theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40),legend.position = "top")+
  scale_fill_brewer(palette = "Accent",direction = -1)+
  facet_grid(.~Cluster)
plot_grid(p1,p2,nrow=2,ncol=1)
dev.off()


#只看肿瘤群体呢？（绝对比例）
TumorClusterRatioInAllCluster<-Cluster_Percet[which(Cluster_Percet$Cluster %in% blast),]
pdf(file = paste(patientID,"_",dim,"_Add_ClusterPercetAndCellCycleChanges_TumorLike_res",res,".pdf",sep=""),width=45,height = 25)
p1<-ggplot(data=TumorClusterRatioInAllCluster, aes(x=SampleID, y=percent, fill=Cluster))+
  geom_bar(stat="identity",position="fill")+
  #	geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
  labs(title=paste(dim,"clusterPercentChanges_res_",res,sep=""))+
  theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40),legend.position = "top")+
#  scale_fill_brewer(palette = "Set3",direction = -1)+
  scale_fill_manual(values = getPalette(colourCount)) 
#  facet_grid(.~Cluster)
plot(p1)
dev.off()

pdf(file = paste(patientID,"_",dim,"_Add_ClusterCellCycleChanges_TumorLike_res",res,".pdf",sep=""),width=45,height = 25)
p1<-ggplot(data=TumorClusterRatioInAllCluster, aes(x=SampleID, y=percent, fill=Phase))+
  geom_bar(stat="identity",position="fill")+
  #	geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
  labs(title=paste(dim,"clusterPercentChanges_res_",res,sep=""))+
  theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40),legend.position = "top")+
#  scale_fill_brewer(palette = "Set3",direction = -1)+
  scale_fill_manual(values = getPalette(colourCount)) 
#  facet_grid(.~Cluster)
plot(p1)
dev.off()


#只看肿瘤群体呢？在PB中的变化
cellInPB <- SampleAlignmentSeurat@meta.data
TumorClusterRatioInAllCluster_PB <- table(cellInPB$stim,cellInPB$orig.ident)
TumorClusterRatioInAllCluster_PB <- prop.table(TumorClusterRatioInAllCluster_PB,margin=1)*100
TumorClusterRatioInAllCluster_PB <- as.data.frame(TumorClusterRatioInAllCluster_PB)
colnames(TumorClusterRatioInAllCluster_PB) <- c("SampleID","Cluster","percent")
TumorClusterRatioInAllCluster_PB <- TumorClusterRatioInAllCluster_PB[which(TumorClusterRatioInAllCluster_PB$Cluster %in% blast),]
TumorClusterRatioInAllCluster_PB$SampleID <- as.numeric(TumorClusterRatioInAllCluster_PB$SampleID)
#TumorClusterRatioInAllCluster_PB = ddply(TumorClusterRatioInAllCluster_PB,"SampleID",transform,percent = percent/sum(percent) * 100)
colourCount = length(unique(TumorClusterRatioInAllCluster_PB$Cluster))
getPalette = colorRampPalette(colors = brewer.pal(9, "Set1"))

pdf(file = paste(patientID,"_",dim,"TumorClusterRatioInAllCluster_PB_res",res,".pdf",sep=""),width=30,height = 25)
p1<-ggplot(data=TumorClusterRatioInAllCluster_PB, aes(x=SampleID, y=percent, fill=Cluster,order = desc(Cluster)))+
  geom_area(stat="identity",alpha=0.6) + #colour = "black",size=0.5,
  scale_fill_manual(values = getPalette(colourCount)) +
  labs(title=paste(patientID,"_",dim,"TumorClusterRatioInAllCluster_PB_res_",res,sep=""))+
  scale_x_continuous(labels=unique(cellInPB$stim)) +
  theme(axis.text.x = element_text(size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40))
plot(p1)
dev.off()

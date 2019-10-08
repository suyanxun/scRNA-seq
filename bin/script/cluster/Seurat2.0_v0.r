#!/bin/sh
#PBS -q core24
#PBS -l mem=30gb,walltime=999:00:00,nodes=1:ppn=10
#HSCHED -s scRNAseq+cellranger+human
#PPN limit 10
#PBS -N merge3sample_10cc
#PBS -o merge3sample_10cc.out
#PBS -e merge3sample_10cc.err
#===============================================================================================================
export  R_LIBS=/pnas/wangqf_group/zhangchen/hanxueshuai/R-3.5.0/library/tibble/libs:$R_LIBS

cd /asnas/wangqf_group/huangyue/seurat/sh
R=/pnas/wangqf_group/zhangchen/hanxueshuai/R-3.5.0/bin
$R/Rscript /home/huangyue/seurat/sh/test11_18.R


#################################################################################################
##                                         part 1                 40min                        ##
#################################################################################################
#remove the R0D and do the sample alignment
rm(list = ls())
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat,lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R/OldVersion"))
library(knitr)
library(dplyr)
library(magrittr)

#设置工作路径
setwd("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/2019Cell/2.0")

#命名
rawSampleName<-c("P101SR00-PB","P105SR00-PB","P106SR00-PB")
sampleName<-c("P101SR00","P105SR00","P106SR00")

sampleList_seurat=c()
sampleList_10x=c()
for (i in 1:length(sampleName)) {
	a<-paste("Sample",i,sep="")
	sampleList_seurat<-c(sampleList_seurat,a)
	b<-paste("Sample",i,"_10X",sep="")
	sampleList_10x<-c(sampleList_10x,b)
}
sampleList_seurat
sampleList_10x

#循环多个样本预处理
for  (j in 1:length(sampleName)) {
    Path<-paste("/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/",rawSampleName[j],"_new_cellno6000/outs/filtered_feature_bc_matrix",sep='')  #导入数据
    list.files(Path)
    Sample_10x<-Read10X(data.dir = Path)
    Sample_10x <- as.matrix(Sample_10x)
    a=c()
    for (i in 1:(length(colnames(Sample_10x)))) {
        a=c(a,paste(sampleName[j],colnames(Sample_10x)[i],sep = "_"))
    }
    colnames(Sample_10x)  <- a
    
    ##该函数还具有数据过滤的功能,基因至少在min.cells个细胞中表达，每个细胞中至少表达min.genes个基因,可以根据具体情况决定这些阈值
    Sample<- CreateSeuratObject(raw.data= Sample_10x,min.cells=3,min.genes=10,project = sampleName[j])
    Sample@meta.data$stim <- sampleName[j]
    ##画图确定阈值
    mito.genes<- grep(pattern = "^MT-",x = rownames(Sample@data),value=T)   #提取线粒体基因列表
    percent.mito <- Matrix::colSums(Sample@raw.data[mito.genes,])/Matrix::colSums(Sample@raw.data)  #计算线粒体基因含量
    Sample <- AddMetaData(object = Sample, metadata = percent.mito, col.name = "percent.mito")  #将线粒体基因含量添加到meta.data>中
    #pdf(file = paste(sampleName[j],"_GeneUMIMito.pdf",sep=''),width=14,height = 8)
    #VlnPlot(object = Sample, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)  #绘制小提琴图
    #dev.off()
    #pdf(file = paste(sampleName[j],"_relationships.pdf",sep=''),width=14,height = 8)    #分析2个指标之间的相互关系
    #par(mfrow = c(1, 2))
    #GenePlot(object = Sample, gene1 = "nUMI", gene2 = "percent.mito")
    #GenePlot(object = Sample, gene1 = "nUMI", gene2 = "nGene")
    #dev.off()
	
    #滤掉小于low.thresholds或大于high.thresholds的数据，这个参数值对应subset.names的参数
    Sample <- FilterCells(object = Sample, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.15))
    #数据标准化，即表达量的计算
	Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize",scale.factor = 10000)
    #计算高度变化的基因（并不是所有基因都有有效的信息，寻找高度变化的基因既可以节省计算资源，又能概括样本信息）
    #pdf(file = paste(sampleName[j],"_FindVariableGenes.pdf",sep=''),width=14,height = 8)
    Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.5)
    #dev.off()
    #cell cycle regression
    cc.genes <- readLines(con = "/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/regev_lab_cell_cycle_genes.txt")
    s.genes <- cc.genes[1:43]
    g2m.genes <- cc.genes[44:97]
    Sample <- CellCycleScoring(object = Sample, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
    Sample@meta.data$CC.Difference <- Sample@meta.data$S.Score - Sample@meta.data$G2M.Score
    #进一步降低系统误差
    Sample <- ScaleData(object = Sample, vars.to.regress = c("nUMI", "percent.mito","CC.Difference"))
    assign(sampleList_seurat[j],Sample)
    rm(Sample_10x)
}


##==============================sample alignment================================
# Gene selection for input to CCA
hvg.union <- c()
for (i in 1:length(sampleList_seurat)) {
  hvg.union <- c(hvg.union, head(rownames(get(sampleList_seurat[i])@hvg.info), 1000))
}

hvg.union <- names(which(table(hvg.union) > 1))
for (i in 1:length(sampleList_seurat)) {
  hvg.union <- hvg.union[hvg.union %in% rownames(get(sampleList_seurat[i])@scale.data)]
}

#run a canonical correlation analysis to identify common sources of variation
SampleAlignment.list <- list(get(sampleList_seurat[1]),get(sampleList_seurat[2]),get(sampleList_seurat[3]))
#时间较长（20cc约10min）
SampleAlignment <- RunMultiCCA(object.list=SampleAlignment.list,genes.use = hvg.union,num.ccs = 20)
#SampleAlignment <- RunMultiCCA(object.list=SampleAlignment.list,genes.use = hvg.union,num.ccs = 20)
#SampleAlignment <- RunMultiCCA(object.list=SampleAlignment.list,genes.use = hvg.union,num.ccs = 20)
#SampleAlignment <- RunMultiCCA(object.list=SampleAlignment.list,genes.use = hvg.union,num.ccs = 20)

# CC Selection类似于PCA分析的PC选择，执行RunCCA后也需要选择合适的CC个数进行后续分析
#pdf(file = paste("SampleAlignment_CCADetermination.pdf"),width=14,height = 8)
#MetageneBicorPlot(SampleAlignment, grouping.var = "stim", dims.eval = 1:20, display.progress = FALSE)
#dev.off()
#pdf(file = paste("SampleAlignment_CCADetermination_heatmap.pdf"),width=14,height = 30)
#DimHeatmap(object = SampleAlignment, reduction.type = "cca", cells.use = 500, dim.use = 1:35,do.balanced = TRUE)
#dev.off()

#saveRDS(SampleAlignment, file = "./SampleAlignment.rds")

#################################################################################################
##                                         part 2                       70min                  ##
#################################################################################################
rm(list = ls())
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat,lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R/OldVersion"))
library(knitr)
library(dplyr)
library(magrittr)

#设置工作路径
setwd("/asnas/wangqf_group/huangyue/seurat/20cc_newout")
sampleName<-c("P101SR00","P105SR00","P106SR00")
res=c(0.6,1,1.4,1.8,2.3,2.7,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4) ##resolution
dim <-15
#flowCytometryFeatures<-c("PTPRC","CD19","MME","MS4A1","CD27","CR2","CD38","IL3RA","ITGAX","HLA-DRA","HLA-DRB1","CD3D","CD3E","CD3G","CD14","FCGR3A","THBD","CD1C","CD34","KIT","CD33","ANPEP","ITGAM","FCGR1A","NCAM1","CD7","CD8A","CD8B","CCR7","CD69","CD4","ITGAE")

#align the CCA subspaces
SampleAlignment <- AlignSubspace(SampleAlignment, reduction.type = "cca", grouping.var = "stim",dims.align = 1:dim) #20cc需要25min
SampleAlignment <- RunTSNE(object = SampleAlignment, reduction.use = "cca.aligned", dims.use = 1:dim,do.fast = TRUE,perplexity = 100)
SampleAlignment <- FindClusters(object = SampleAlignment, reduction.type = "cca.aligned", dims.use = 1:dim,save.SNN = TRUE,resolution = res) #20cc需要40min

saveRDS(SampleAlignment, file = "./SampleAlignment_merge3_20.rds")


#################################################################################################
##                                         part 3                                              ##
#################################################################################################
rm(list = ls())
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat)
library(knitr)
library(dplyr)
library(magrittr)

dim <-22
res=c(0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.2,3.5,4) ##resolution
setwd("/home/huangyue/seurat/out")
SampleAlignment<-readRDS(file = "./SampleAlignment_1_15.rds")

####整个res一起画
p<-c()
for (i in 1:(dim(SampleAlignment@meta.data)[2]-9)) {
    a<-paste("p",i,sep="")
    p<-c(p,a)
}
for (i  in c(10:(dim(SampleAlignment@meta.data)[2]))) {
    temp <- TSNEPlot(object = SampleAlignment, group.by = colnames(SampleAlignment@meta.data[i]),plot.order=c(levels(factor(SampleAlignment@meta.data[i],levels=c((dim(unique(SampleAlignment@meta.data[i]))[1]-1):0)))),do.return = T ,pt.size = 0.5,do.label = TRUE,label.size=8, plot.title=colnames(SampleAlignment@meta.data[i]))
    assign(p[i-9],temp)
}

png(file = paste(dim,"TSNEPlot_differentRes.png",sep=""),width=4500,height = 3000)
plot_grid(get(p[1]),get(p[2]),get(p[3]),get(p[4]),get(p[5]),get(p[6]),get(p[7]),get(p[8]), get(p[9]),get(p[10]),get(p[11]),get(p[12]),get(p[13]),get(p[14]),get(p[15]))
dev.off()
##单个res
#Visualization
#png(file = paste(dim,"Samplemerge_TSNEPlot_res",res,".png",sep=""),width=1350,height = 800)
#p1 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P101SR00","P105SR00","P106SR00"),do.return = TRUE, pt.size = 0.2) #different sample
#p2 <- TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 0.5,do.label = TRUE,label.size=8) #different cluster
#plot_grid(p1, p2)
#dev.off()


features<-c("MTOR","PTEN","NCAM1","KIT","ITGAM","MPL","PTPRC","FLT3","CD7","MME","CD34","CD38","IL3RA","PROM1","CD79A","MS4A1","IGHD","CD3D","CD3E","CD3G","IL7R","CD4","CD8A","CD8B","GZMB","PRF1","TNFRSF18","CCR10","ID3","FCGR3A","NKG7","GNLY","KLRC1","CD14","LYZ","S100A9","CSF3R","MS4A7","CD33","CCL2","LYN","CSF1R","IFITM1","IFITM2","IFITM3","S100A8","CD1C","ITGAX","CLEC4C","GP9","ITGA2B","PF4","PPBP","CX3CR1","TMEM40","HBA1","HBA2","HBB","GYPA","AZU1","ELANE","CEACAM8","LY6G5B","LY6G5C")


ncol=7
nrow= ceiling(length(features)/ncol)

png(file = paste(dim,"cluster_marker_3sample_15cc.png",sep=""),width=675*ncol,height = 800*nrow,pointsize=10)
FeaturePlot(object = SampleAlignment, features.plot = features, cols.use = c("#CCCCCC", "#FF0000"), reduction.use = "tsne",pt.size=1,nCol=ncol)
dev.off()
####可以画特定res的图
png(file = paste(dim,"VlnPlot_marker_",colnames(SampleAlignment@meta.data[22]),".png",sep=""),width=3360,height =320*nrow ,pointsize=10)
VlnPlot(object = SampleAlignment, point.size.use = 0,group.by = colnames(SampleAlignment@meta.data[22]),features.plot =features)
dev.off()
####可以画特定样本的图
#SampleAlignment101<- SubsetData(object=SampleAlignment, cells.use =celluse101)
#SampleAlignment105<- SubsetData(object=SampleAlignment, cells.use =celluse105)
#SampleAlignment106<- SubsetData(object=SampleAlignment, cells.use =celluse106)
#SampleAlignment101<-SetIdent(SampleAlignment101,cells.use = NULL, ident.use =SampleAlignment101@meta.data$cellident)
#png(file = paste("VlnPlot_marker_101.png"),width=3360,height =320*nrow,pointsize=10)
#VlnPlot(object = SampleAlignment101, point.size.use = 0.5,features.plot =features)
#dev.off()


#find the cluster marker in every cluster
ClusterMarker<-FindAllMarkers(object = SampleAlignment,only.pos = T,min.pct = 0.25,thresh.use=0.25)  #15cc大约100min
write.table(ClusterMarker,file = paste(dim,"ClusterMarker.txt",sep=""),quote = FALSE,row.names = TRUE,sep = "\t")


ClusterMarker101<-FindAllMarkers(object = SampleAlignment101,only.pos = T,min.pct = 0.25,thresh.use=0.25)  #15cc大约100min
write.table(ClusterMarker101,file = paste("ClusterMarker101.txt",sep=""),quote = FALSE,row.names = TRUE,sep = "\t")
ClusterMarker105<-FindAllMarkers(object = SampleAlignment105,only.pos = T,min.pct = 0.25,thresh.use=0.25)  #15cc大约100min
write.table(ClusterMarker105,file = paste("ClusterMarker105.txt",sep=""),quote = FALSE,row.names = TRUE,sep = "\t")


png(file = paste(dim,"heatmap_",res,".png",sep=""),width=3000,height = 10000)
top100 <- ClusterMarker %>% group_by(cluster) %>% top_n(100, avg_logFC)
p<-DoHeatmap(object = SampleAlignment, genes.use = top100$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
p$data
#write.table(top100,paste(dim,"ClusterMarkerTop100_",res,".txt",sep=""),sep="\t")

saveRDS(SampleAlignment, file = "./SampleAlignment_2_15.rds")


#################################################################################################
##                                         part 4            new.ident                         ##
#################################################################################################
rm(list = ls())
library(ggplot2)
library(cowplot)
library(Matrix)
library(Seurat,lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R/OldVersion"))
library(knitr)
library(dplyr)
library(magrittr)

setwd("/asnas/wangqf_group/huangyue/seurat/out")
SampleAlignment<-readRDS(file = "./SampleAlignment_3_15.rds")

res=3.2 ##resolution
dim <- 15

#ident 15cc,3.2res
new.ident <-c("CD14mono_1","stem_1","Granu","CD8+ T","CD8+ T","plasma","stem_3","stem_1","CD4+T","Progenitor","CD14mono_1","CD14mono_1","CD16mono","stem_2","CD4+T","CD4+T","CD4+T","CD14mono_2","cDC","CD14mono_1","CD4+T","stem_2","CD4+T","NK","stem_4","B","Eryth","CD4+T","Granu","CD4+T","MK","stem_1","B","stem_1","stem_3","CD4+T")


"Stem_3","CD14","CD14","CD8+ T","Granu","Granu","B cells","CD14","Stem_3","CD4+ T","CD8+ T","CD16","CD16","DC","Stem_1","Proge_1","Stem_2","CD4+ T","CD4+ T","CD4+ T","CD4+ T","unkown","NK","Proge_2","CD4+ T","Stem_3","Eryth","CD16","CD4+ T","Granu","CD8+ T","Stem_3","MK","Granu","Proge_2","CD4+ T","Stem_3")
#ident_new
new.ident <- c("Stem_1","CD14","Granu","Granu","CD8+ T","B cells","CD14","Granu","CD4+ T","Proge_2","CD8+ T","Proge_2","DC","CD4+ T","CD16","Stem_2","CD14","CD4+ T","Stem_3","CD4+ T","CD4+ T","CD14","CD4+ T","NK","Stem_3","Eryth","Proge_1","CD4+ T","Granu","CD8+ T","CD4+ T","MK","CD8+ T","Stem_3","Proge_1")

for (i in 0:(length(new.ident)-1)) {
    SampleAlignment<- RenameIdent(object =SampleAlignment,old.ident.name = i,new.ident.name = new.ident[i + 1])
}

#画出样本来源图
png(file = paste(dim,"Samplemerge_TSNEPlot.png",sep=""),width=5400,height = 1600)
p1 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P101SR00","P105SR00","P106SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#FF450080","#00CD6680","#00B2EE80")) #different sample
p2 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P106SR00","P105SR00","P101SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#FF450080"))
p3 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P105SR00","P101SR00","P106SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#00CD6680"))
p4 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P101SR00","P105SR00","P106SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#00B2EE80"))
plot_grid(p1, p2,p3,p4,ncol=4)
dev.off()

#画出细胞周期来源图
png(file = paste(dim,"Samplemerge_TSNEPlot.png",sep=""),width=5400,height = 1600)
p1 <- TSNEPlot(object = SampleAlignment, group.by = "Phase",  plot.order=c("G2M","S","G1"),do.return = TRUE, pt.size = 3,colors.use=c("#FF450080","#00B2EE80","#00CD6680")) #different sample
p2 <- TSNEPlot(object = SampleAlignment, group.by = "Phase",  plot.order=c("G1","G2M","S"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#FF450080"))
p3 <- TSNEPlot(object = SampleAlignment, group.by = "Phase",  plot.order=c("G2M","G1","S"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#00CD6680"))
p4 <- TSNEPlot(object = SampleAlignment, group.by = "Phase",  plot.order=c("S","G1","G2M"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#00B2EE80"))
plot_grid(p1, p2,p3,p4,ncol=4)
dev.off()

#画出分群图
png(file = paste(dim,"Samplemerge_TSNEPlot.png",sep=""),width=5400,height = 1600)
p1<-TSNEPlot(object = SampleAlignment,do.return = T ,pt.size = 3,do.label = TRUE,label.size=15,plot.title="merge")
p2<-TSNEPlot(object = SampleAlignment,do.return = T ,pt.size = 3,do.label = TRUE,label.size=15,cells.use=celluse101,plot.title="P101")
p3<-TSNEPlot(object = SampleAlignment,do.return = T ,pt.size = 3,do.label = TRUE,label.size=15,cells.use=celluse105,plot.title="P105")
p4<-TSNEPlot(object = SampleAlignment,do.return = T ,pt.size = 3,do.label = TRUE,label.size=15,cells.use=celluse106,plot.title="P106")
plot_grid(p1, p2,p3,p4,ncol=4)
dev.off()


#show the ratio changes of each cluster
#=================================calculate the cell cycle============================
#cc.genes <- readLines(con = "/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/regev_lab_cell_cycle_genes.txt")
#s.genes <- cc.genes[1:43]
#g2m.genes <- cc.genes[44:97]
#calculate the cell cycle score
#SampleAlignment<- CellCycleScoring(object = SampleAlignment, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

#table(SampleAlignment@meta.data$stim,SampleAlignment@meta.data$Phase) #check the result
#SampleAlignment@meta.data$CC.Difference <- SampleAlignment@meta.data$S.Score - SampleAlignment@meta.data$G2M
#SampleAlignment<- ScaleData(object =SampleAlignment, vars.to.regress = c("CC.Difference")


#Freq<-table(SampleAlignment@meta.data$stim,SampleAlignment@ident,SampleAlignment@meta.data$Phase)
#write.table(Freq,file= paste(dim,"clusterFreq_res",res,".txt",sep=""),quote = FALSE,row.names = FALSE,sep = "\t")
#Freq<-read.table(file=paste(dim,"clusterFreq_res",res,".txt",sep=""),header=T,sep="\t")
#colnames(Freq)<-c("SampleID","Cluster","Phase","Count")
#Cluster_Percet = ddply(Freq,"SampleID",transform,percent = Count/sum(Count) * 100)
#Cluster_Percet$SampleID<-factor(Cluster_Percet$SampleID,levels=c("P101SR00","P105SR00","P106SR00"))
#Cluster_Percet$Phase<-factor(Cluster_Percet$Phase,levels=c("G2M","S","G1"))

#pdf(file = paste(dim,"ClusterPercetAndCellCycleChanges_res",res,".pdf",sep=""),width=90,height = 25)
#p1<-ggplot(data=Cluster_Percet, aes(x=SampleID, y=percent, fill=SampleID))+geom_bar(stat="identity")+
#    geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
#    labs(title=paste(dim,"clusterPercentChanges_res_",res,sep=""))+
#    theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 30),text = element_text(size = 40),legend.position = "top")+
#    scale_fill_brewer(palette = "Set3",direction = -1)+
#   facet_grid(.~Cluster)
#p2<-ggplot(data=Cluster_Percet, aes(x=SampleID, y=percent, fill=Phase))+
#   geom_bar(stat="identity",position="fill")+
#   geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
#   labs(title=paste(dim,"clusterCellCycleChanges_res_",res,sep=""))+
#    theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 30),text = element_text(size = 40),legend.position = "top")+
#    scale_fill_brewer(palette = "Accent",direction = -1)+
#   facet_grid(.~Cluster)
#plot_grid(p1,p2,nrow=2,ncol=1)
#dev.off()

#saveRDS(SampleAlignment, file = "SampleAlignment3.rds")


############################################################################################################################
#=================================Visualization=================================
#Visualization with new ident
#UMICount<-read.table(file="/asnas/wangqf_group/jiangsht/SingleCellRNAseq/Result/IntegratingMultiSample/Test/Test4CountMutation/ExtractCBFromScRNA/PysamExtractReads/P101SRMutationDetection_KIT.txt",header=TRUE)
#test<-subset(UMICount,select=c("finalCellMutaionStatus"))
#colnames(test)<-paste(colnames(test),"EZH2",sep="_")
#SampleAlignmentSeurat<-AddMetaData(object = SampleAlignmentSeurat, metadata = test)

#################################################################################################
##                                         part 5            metadata                         ##
#################################################################################################
#metadata
cellident<-SampleAlignment@cell.names
mutationStatus<-rep("N",times=length(cellident))
finalMutation<-data.frame(cellident,mutationStatus)
rownames(finalMutation)<-finalMutation[,1]
finalMutation$mutationStatus<-as.character(finalMutation$mutationStatus)

#table101_EZH2
p101_EZH2<-read.table("/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/P101SRMutationDetection_EZH2.txt",head=TRUE)
test <-subset(p101_EZH2,select=c("patientIDCB","finalCellMutaionStatus"))
test$finalCellMutaionStatus<-paste(test$finalCellMutaionStatus,"EZH2",sep="_")

overlapCell<-which(rownames(finalMutation) %in% rownames(test))
for (i in overlapCell) {
    finalMutation[i,"mutationStatus"]<-test[rownames(finalMutation)[i],"finalCellMutaionStatus"]
}

rownames(subset(finalMutation,mutationStatus=="Mut_EZH2"))
table(finalMutation$mutationStatus)

#table101_KIT
p101_KIT<-read.table("/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/P101SRMutationDetection_KIT.txt",head=TRUE)
test2 <-subset(p101_KIT,select=c("finalCellMutaionStatus"))
test2$finalCellMutaionStatus<-paste(test2$finalCellMutaionStatus,"KIT",sep="_")
#寻找KIT的（查看是否重叠），填写好替换
#finalMutation2<-finalMutation
#for (i in 1:dim(finalMutation2)[1]) {
#   finalMutation2[i,"mutationStatus"]<-test2[rownames(finalMutation2)[i],"finalCellMutaionStatus"]
#}
#rownames(subset(finalMutation2,mutationStatus=="Mut_KIT"))
"P101SR00_TACGGTAGTGACTACT" "P101SR00_TCTGAGAAGTACGCCC"
"P101SR00_CTAAGACGTCTTCGTC" "P101SR00_GAATGAAAGCACGCCT"
finalMutation["P101SR00_TACGGTAGTGACTACT","mutationStatus"]<-test2["P101SR00_TACGGTAGTGACTACT","finalCellMutaionStatus"]
finalMutation["P101SR00_TCTGAGAAGTACGCCC","mutationStatus"]<-test2["P101SR00_TCTGAGAAGTACGCCC","finalCellMutaionStatus"]
finalMutation["P101SR00_CTAAGACGTCTTCGTC","mutationStatus"]<-test2["P101SR00_CTAAGACGTCTTCGTC","finalCellMutaionStatus"]
finalMutation["P101SR00_GAATGAAAGCACGCCT","mutationStatus"]<-test2["P101SR00_GAATGAAAGCACGCCT","finalCellMutaionStatus"]

rownames(subset(finalMutation,mutationStatus=="Mut_KIT"))
table(finalMutation$mutationStatus)

#table105_PTPN11
p105_PTPN11<-read.table("/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/P105SRMutationDetection.txt",head=TRUE)
test3<-subset(p105_PTPN11,patientID=="P105SR00-PB")
test3 <-subset(test3,select=c("patientIDCB","finalCellMutaionStatus"))
test3$finalCellMutaionStatus<-paste(test3$finalCellMutaionStatus,"PTPN11",sep="_")
###修改读取的数据patientIDCB列
test3$patientIDCB<-as.character(test3$patientIDCB)
for (i in c(1:9131)){
    rownames(test3)[i]<-paste(substr(rownames(test3)[i],1,8),substr(rownames(test3)[i],12,28),sep="")
}
head(test3)

#查找105各点
#finalMutation3<-finalMutation
#overlapCell<-which(rownames(finalMutation3) %in% rownames(test3))
#for (i in overlapCell) {
#    finalMutation3[i,"mutationStatus"]<-test3[rownames(finalMutation3)[i],"finalCellMutaionStatus"]
#}

#rownames(subset(finalMutation3,mutationStatus=="WT_PTPN11"))
total<-c("P105SR00_ACCTTTATCAGCACAT","P105SR00_CATTCGCCACTTCTGC","P105SR00_TCTGGAACACATTTCT","P105SR00_TGAGGGAGTTTGTTTC","P105SR00_TTATGCTGTGGTACAG","P105SR00_AACTCAGTCTCCAGGG","P105SR00_ACGTCAAGTCGAAAGC","P105SR00_AGTCTTTGTGATGTGG","P105SR00_ATAAGAGTCTCCTATA","P105SR00_CCCTCCTAGTAAGTAC","P105SR00_CGGAGCTGTTCACGGC","P105SR00_CGTCAGGGTCAGAATA","P105SR00_CGTGTCTCACAGCGTC","P105SR00_CTGCGGAAGTACTTGC","P105SR00_GACCTGGGTTGTCTTT","P105SR00_TACCTTAGTCATATGC","P105SR00_TACTTGTCACACGCTG")
for (i in total){
    finalMutation[i,"mutationStatus"]<-test3[i,"finalCellMutaionStatus"]
}
table(finalMutation$mutationStatus)

SampleAlignment<-AddMetaData(object = SampleAlignment, metadata = finalMutation)


SampleAlignment_merge_101<-readRDS(file ="/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/FiveSampleAlignment_13.rds")
SampleAlignment_merge_105<-readRDS(file ="/asnas/wangqf_group/jiangsht/Share/MutationInscRNA/P105SR_10SampleAlignmentSeurat.rds")
SampleAlignment_merge_106<-readRDS(file ="/asnas/wangqf_group/huangyue/seurat/15cc_newout/SampleAlignment_only106.rds")
###点亮推断初诊肿瘤细胞105
cellident<-SampleAlignment@cell.names
origident<-rep("N",times=length(cellident))
finalorigident<-data.frame(cellident,origident)
finalorigident$origident<-as.character(finalorigident$origident)
rownames(finalorigident)<-finalorigident[,1]

test105 <-subset(SampleAlignment_merge_105@meta.data,select=c("orig.ident","finalCellMutaionStatus_PTPN11","stim"))
colnames(test105)[1]<-paste("origident")
test105<-test105[test105$stim=="P105SR00-PB",]
test105$finalCellMutaionStatus_PTPN11<-paste(test105$finalCellMutaionStatus_PTPN11,"PTPN11_check",sep="_")
rownames(test105)<-as.character(rownames(test105))
#length()
for (i in c(1:9131)){
    rownames(test105)[i]<-paste(substr(rownames(test105)[i],1,8),substr(rownames(test105)[i],12,28),sep="")
}
head(test105)
test105$origident<-as.character(test105$origident)
overlapCell<-which(rownames(finalorigident) %in% rownames(test105))
for (i in overlapCell) {
    finalorigident[i,"origident"]<-test105[rownames(finalorigident)[i],"origident"]
}


table(finalorigident$origident)

###点亮推断初诊肿瘤细胞101

origident101<-SampleAlignment_merge_101@cell.names
test101 <-subset(SampleAlignment_merge_101@meta.data,select=c("orig.ident","finalCellMutaionStatus_EZH2"))
colnames(test101)[1]<-paste("origident101")
test101$finalCellMutaionStatus_EZH2<-paste(test101$finalCellMutaionStatus_EZH2,"EZH2_check",sep="_")
rownames(test101)<-as.character(rownames(test101))
test101$origident101<-as.character(test101$origident101)
overlapCell<-which(rownames(finalorigident) %in% rownames(test101))
for (i in overlapCell) {
    finalorigident[i,"origident"]<-test101[rownames(finalorigident)[i],"origident101"]
}

###点亮推断初诊肿瘤细胞106

origident106<-SampleAlignment_merge_106@cell.names
test106 <-subset(SampleAlignment_merge_106@meta.data,select=c("cellident","res.4"))
colnames(test106)[1]<-paste("origident106")
test106$res.4<-paste("N")
rownames(test106)<-as.character(rownames(test106))
test106$origident106<-as.character(test106$origident106)
##修改一致的rownames
for (i in c(1:8346)){
    rownames(test106)[i]<-paste(substr(rownames(test106)[i],1,8),substr(rownames(test106)[i],12,28),sep="")
}

overlapCell<-which(rownames(finalorigident) %in% rownames(test106))
for (i in overlapCell) {
    finalorigident[i,"origident"]<-test106[rownames(finalorigident)[i],"origident106"]
}

SampleAlignment<-AddMetaData(object = SampleAlignment, metadata = finalorigident)
table(finalorigident$origident)


p3 <- TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_PTPN11_check","Mut_PTPN11_check","N","NoCoverage_PTPN11_check"),do.return = TRUE, pt.size = 5,colors.use=c("#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"))


celluse<-grep("P105SR00",rownames(SampleAlignment@meta.data),value=T)
png(file = paste(dim,"TSNEPlot.png",sep=""),width=1350,height = 1600)
TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("stemCell_1","stemCell_2","Proge_1","unknown","B","CD14Mono","CD16Mono_1","CD16Mono_2","CD4 Treg","CD4 memo T","CD4T","CD8 memo T","CD8T","Eryth","Granu","MK","N","NK","NKT","Proge_2","Proge_3","cDC","pDC"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse,do.return = T ,pt.size = 3,do.label = TRUE,label.size=15)
dev.off()

SampleAlignment<-AddMetaData(object = SampleAlignment, metadata = finalorigident)




###cellident添加至metadata
ident<-as.character(SampleAlignment@ident)
SampleAlignment@meta.data$cellident<-as.character(SampleAlignment@meta.data$cellident)
for (i in c(1:length(SampleAlignment@ident))){
    SampleAlignment@meta.data$cellident[i]<-ident[i]
}
head(SampleAlignment@meta.data$cellident)

png(file = paste(dim,"TSNEPlot_merge.png",sep=""),width=5400,height = 1600)
#pdf(file = paste(dim,"Samplemerge_cellCycle_TSNEPlot_res",res,".pdf",sep=""),width=25,height = 8)
p1 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=names(sort(table(SampleAlignment@meta.data$stim))),do.return = TRUE, pt.size = 3) #different sample
p2 <- TSNEPlot(object = SampleAlignment, do.return = TRUE ,pt.size = 3,do.label = TRUE,label.size=10,plot.title ="TSNEPlot") #different cluster
p3 <- TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("Mut_EZH2","Mut_KIT","Mut_PTPN11","N","NoCoverage_EZH2","WT_EZH2","WT_KIT","WT_PTPN11","further2confirm_EZH2"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#3A5FCD","#CCCCCC","#CCCCCC","#CCCCCC","#EE0000","#CCCCCC"))

p1 <- TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_EZH2","Mut_EZH2","Mut_KIT","Mut_PTPN11","N","NoCoverage_EZH2","WT_KIT","WT_PTPN11","further2confirm_EZH2"),do.return = TRUE, pt.size = 5,colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"))
p2 <- TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_KIT","Mut_KIT","WT_EZH2","Mut_EZH2","Mut_PTPN11","N","NoCoverage_EZH2","WT_PTPN11","further2confirm_EZH2"),do.return = TRUE, pt.size = 5,colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"))
p3 <- TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_PTPN11","Mut_PTPN11","WT_KIT","Mut_KIT","WT_EZH2","Mut_EZH2","N","NoCoverage_EZH2","further2confirm_EZH2"),do.return = TRUE, pt.size = 5,colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"))
#p4 <- TSNEPlot(object = SampleAlignment, group.by = "mutationStatus",plot.order=c("WT"," Mut_KIT","NoCoverage"),do.return = TRUE, pt.size = 2,colors.use=c("#CCCCCC","#CC0033","#0000FF"))
#p5 <- TSNEPlot(object = SampleAlignment, group.by = "Phase",do.return = TRUE, pt.size = 1)
plot_grid(p1,p2,p3,ncol=3)
dev.off()

##merge_105.png
png(file = paste(dim,"merge_105_unknown.png",sep=""),width=5400,height = 1600)
p1 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P101SR00","P105SR00","P106SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#FF450080","#00CD6680","#00B2EE80"))
p2 <- TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 3,do.label = TRUE,label.size=10)
celluse<-grep("P105SR00",rownames(SampleAlignment@meta.data),value=T)
p3 <-TSNEPlot(object = SampleAlignment,group.by = "origident",cells.use=celluse,do.return = T ,pt.size = 3,do.label = TRUE,label.size=10)
p4 <-TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("stemCell_2","stemCell_1","Proge_1","unknown","B","CD14Mono","CD16Mono_1","CD16Mono_2","CD4 Treg","CD4 memo T","CD4T","CD8 memo T","CD8T","Eryth","Granu","MK","N","NK","NKT","Proge_2","Proge_3","cDC","pDC"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#FFC125"),cells.use=celluse,do.return = T ,pt.size = 3,do.label = TRUE,label.size=10)

plot_grid(p1,p2,p3,p4,ncol=4)
plot_grid(p1,p2,p3,p4,ncol=4)
dev.off()

saveRDS(SampleAlignment, file = "./SampleAlignment_3_15.rds")


###计算临床blast
#1.判断突变在哪些群体中，定义为肿瘤细胞群体
#2.查看各肿瘤群体的数目
table(SampleAlignment@ident)
a<-subset(SampleAlignment@ident,SampleAlignment@meta.data$stim=="P105SR00")

####11111聚类查看群体间相似性
newname<-paste(SampleAlignment@meta.data$cellident,SampleAlignment@meta.data$stim,sep="_")
newname<-data.frame(newname)
rownames(newname)<-SampleAlignment@cell.names
newname$newname<-as.character(newname$newname)
SampleAlignment<-AddMetaData(object = SampleAlignment, metadata = newname)
SampleAlignment<-SetIdent(SampleAlignment,cells.use = NULL, ident.use =SampleAlignment@meta.data$newname)

library(pheatmap)
a<-AverageExpression(object =SampleAlignment)
png(file = paste("pheatmap_2.png"),width=2000,2000)
pheatmap(a,scal="row",show_rownames=F,fontsize=25)
dev.off()

####22222聚类查看群体间相似性-指定基因,指定ident（通过SetIdent修改）
ClusterMarker<-FindAllMarkers(object = SampleAlignment,only.pos = T,min.pct = 0.25,thresh.use=0.25)  #9:54
write.table(ClusterMarker,file = paste("AllClusterMarker_3sample_15cc.txt",sep=""),quote = FALSE,row.names = TRUE,sep = "\t")
saveRDS(SampleAlignment, file = "./ClusterMarker_15cc.rds")
ClusterMarker<-read.table("/asnas/wangqf_group/huangyue/seurat/15cc_newout/AllClusterMarker_3sample_15cc.txt",head=TRUE,blank.lines.skip=F, sep='\t')
top100 <- ClusterMarker %>% group_by(cluster) %>% top_n(100, avg_logFC)
p<-DoHeatmap(object = SampleAlignment, genes.use = top100$gene, do.plot =F, remove.key = TRUE)
head(p$data)
m<-unique(p$data$gene)
average<-AverageExpression(object =SampleAlignment,genes.use=m)
#matrix=cor(average)
#test101<-subset(average,select=c("B_P101SR00","CD4+T_P101SR00","CD8+ T_P101SR00","CD14mono_1_P101SR00","CD14mono_2_P101SR00","CD16mono_P101SR00","cDC_P101SR00","Eryth_P101SR00","Granu_P101SR00","MK_P101SR00","NK_P101SR00","plasma_P101SR00","Progenitor_P101SR00","stem_1_P101SR00","stem_2_P101SR00","stem_3_P101SR00","stem_4_P101SR00"))
#test105<-subset(average,select=c("B_P105SR00","CD4+T_P105SR00","CD8+ T_P105SR00","CD14mono_1_P105SR00","CD14mono_2_P105SR00","CD16mono_P105SR00","cDC_P105SR00","Eryth_P105SR00","Granu_P105SR00","MK_P105SR00","NK_P105SR00","plasma_P105SR00","Progenitor_P105SR00","stem_1_P105SR00","stem_2_P105SR00","stem_3_P105SR00","stem_4_P105SR00"))
#test106<-subset(average,select=c("B_P106SR00","CD4+T_P106SR00","CD8+ T_P106SR00","CD14mono_1_P106SR00","CD14mono_2_P106SR00","CD16mono_P106SR00","cDC_P106SR00","Eryth_P106SR00","Granu_P106SR00","MK_P106SR00","NK_P106SR00","plasma_P106SR00","Progenitor_P106SR00","stem_1_P106SR00","stem_2_P106SR00","stem_3_P106SR00","stem_4_P106SR00"))
png(file = paste("pheatmap101.png"),width=2000,2000)
annotation_col = data.frame(CellType = factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Tumor","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Tumor","Tumor","Normal","Tumor","Tumor","Normal","Tumor","Tumor","Normal","Tumor","Tumor","Normal")))
rownames(annotation_col)=colnames(average)
pheatmap(average,scal="row",cluster_rows = F,show_rownames=F,fontsize=15,annotation_col=annotation_col)
dev.off()


###指定D0整个群体颜色  101突变
celluse101<-grep("P101SR00",rownames(SampleAlignment@meta.data),value=T)
png(file = paste(dim,"101_EZH2.png",sep=""),width=2700,height = 1600)
p1<-TSNEPlot(object = SampleAlignment,group.by = "cellident",plot.order=c("Stem_1","Stem_2","Stem_3","Granu","Proge_2","Proge_1","B cells","CD14","CD16","CD4+ T","CD8+ T","DC","Eryth","MK","NK"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#B452CD","#00B2EE","#008B00","#FF7F00","#EE4000"),do.return = T ,pt.size = 3,do.label = TRUE,label.size=12,cells.use=celluse101,plot.title="P101_tumor")
p2<-TSNEPlot(object =SampleAlignment,group.by="mutationStatus",plot.order=c("WT_EZH2","Mut_EZH2","1","3","NoCoverage_EZH2","further2confirm_EZH2","N"),do.return = TRUE, pt.size = 1,colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#EE0000","#3A5FCD"),cells.use=celluse101,plot.title="P101_EZH2")
plot_grid(p1,p2)
dev.off()

###105对回D0点亮初诊  newnew
celluse105<-grep("P105SR00",rownames(SampleAlignment@meta.data),value=T)
png(file = paste(dim,"105_D0.png",sep=""),width=2700,height = 1600)
p1<-TSNEPlot(object = SampleAlignment,group.by = "cellident",plot.order=c("Stem_1","Stem_2","Stem_3","Granu","Proge_2","Proge_1","B cells","CD14","CD16","CD4+ T","CD8+ T","DC","Eryth","MK","NK"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#B452CD","#00B2EE","#008B00","#FF7F00","#EE4000"),do.return = T ,pt.size = 3,do.label = TRUE,label.size=12,cells.use=celluse105,plot.title="P105")
p2<-TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("stemCell_1","stemCell_2","stemCell_3","pro_1","CD14Mono","cDC","B_1","B_2","CD16Mono_1","CD16Mono_2","CD4 memo T","CD4 Treg","CD4T","CD8 memo T","CD8T","Eryth","Granu_1","Granu_2","MK","N", "NK","NKT","pDC","pro_2","pro_3","stemCell_4"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#00B2EE","#B452CD","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse105,do.return = T ,pt.size = 3,do.label = TRUE,label.size=8,plot.title="P105_tumor")
plot_grid(p1,p2)
dev.off()

###101对回D0点亮初诊  newnew
celluse101<-grep("P101SR00",rownames(SampleAlignment@meta.data),value=T)
png(file = paste(dim,"101_D0.png",sep=""),width=2700,height = 1600)
p1<-TSNEPlot(object = SampleAlignment,group.by = "cellident",plot.order=c("Stem_1","Stem_2","Stem_3","Proge_2","Granu","Proge_1","B cells","CD14","CD16","CD4+ T","CD8+ T","DC","Eryth","MK","NK"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#43CD80","#FF82AB","#00B2EE","#008B00","#FF7F00","#EE4000"),do.return = T ,pt.size = 3,do.label = TRUE,label.size=12,cells.use=celluse101,plot.title="P101")
p2<-TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("CD34_1","CD34_2","CD34_3","GMP_1","GMP_2","Granulocyte","unknown","B","CD14Mono","CD16Mono","CD4T","CD8T","DC","Eryth","MK","NK","N"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#00B2EE","#FF82AB","#43CD80","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse101,do.return = T ,pt.size = 1,do.label = TRUE,label.size=0,plot.title="P101_tumor")
plot_grid(p1,p2)
dev.off()

###106对回D0点亮初诊  newnew
celluse106<-grep("P106SR00",rownames(SampleAlignment@meta.data),value=T)
png(file = paste(dim,"106_D0.png",sep=""),width=2700,height = 1600)
p1<-TSNEPlot(object = SampleAlignment,group.by = "cellident",plot.order=c("Stem_1","Stem_2","Stem_3","Proge_2","Granu","Proge_1","B cells","CD14","CD16","CD4+ T","CD8+ T","DC","Eryth","MK","NK"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#43CD80","#FF82AB","#00B2EE","#008B00","#FF7F00","#EE4000"),do.return = T ,pt.size = 3,do.label = TRUE,label.size=12,cells.use=celluse106,plot.title="P101")
p2<-TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("Stem","Granu_1","Proge_1","DC","CD14","B","CD16","CD4T","CD8T","Eryth","Granu_2","MK","NK","Proge_2","N"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#00B2EE","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse106,do.return = T ,pt.size = 1,do.label = TRUE,label.size=0,plot.title="P106_tumor")


"CD34_1","CD34_2","CD34_3","GMP_1","GMP_2","Granulocyte","unknown","B","CD14Mono","CD16Mono","CD4T","CD8T","DC","Eryth","MK","NK","N"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#00B2EE","#FF82AB","#43CD80","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse101,do.return = T ,pt.size = 3,do.label = TRUE,label.size=8,plot.title="P101_tumor")
plot_grid(p1,p2)
dev.off()


celluse101<-grep("P101SR00",rownames(SampleAlignment@meta.data),value=T)
celluse105<-grep("P105SR00",rownames(SampleAlignment@meta.data),value=T)
celluse106<-grep("P106SR00",rownames(SampleAlignment@meta.data),value=T)

png(file = paste(dim,"total.png",sep=""),width=6750,height = 1600)
p1 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P101SR00","P105SR00","P106SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#FF450080","#00CD6680","#00B2EE80"))
p2<-TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 3,do.label = TRUE,label.size=10)
p3<-TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 3,do.label = TRUE,label.size=10,cells.use=celluse101)
p4<-TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 3,do.label = TRUE,label.size=10,cells.use=celluse105)
p5<-TSNEPlot(object = SampleAlignment, do.return = T ,pt.size = 3,do.label = TRUE,label.size=10,cells.use=celluse106)
plot_grid(p1,p2,p3,p4,p5,ncol=5)
dev.off()

png(file = paste("pt13.png"),width=1350,height = 1600)
 <- TSNEPlot(object = SampleAlignment, group.by = "stim",  plot.order=c("P106SR00","P101SR00","P105SR00"),do.return = TRUE, pt.size = 3,colors.use=c("#CCCCCC","#CCCCCC","#FF450080"))


pdf(file = paste("pt13.png"),width=18,height = 16)
TSNEPlot(object = SampleAlignment,group.by = "newname",plot.order=c("Granu_P101SR00","Progenitor_P101SR00","stem_1_P101SR00","stem_2_P101SR00","stem_3_P101SR00","stem_4_P101SR00","Granu_P105SR00","Progenitor_P105SR00","stem_1_P105SR00","stem_2_P105SR00","stem_3_P105SR00","stem_4_P105SR00","Granu_P106SR00","Progenitor_P106SR00","stem_1_P106SR00","stem_2_P106SR00","stem_3_P106SR00","stem_4_P106SR00","cDC_P105SR00","cDC_P106SR00","B_P101SR00","B_P105SR00","B_P106SR00","CD14mono_1_P101SR00","CD14mono_1_P105SR00","CD14mono_1_P106SR00","CD14mono_2_P101SR00","CD14mono_2_P105SR00","CD14mono_2_P106SR00","CD16mono_P101SR00","CD16mono_P105SR00","CD16mono_P106SR00","CD4+T_P101SR00","CD4+T_P105SR00","CD4+T_P106SR00","CD8+ T_P101SR00","CD8+ T_P105SR00","CD8+ T_P106SR00","Eryth_P101SR00","Eryth_P105SR00","Eryth_P106SR00","MK_P101SR00","MK_P105SR00","MK_P106SR00","NK_P101SR00","NK_P105SR00","NK_P106SR00","cDC_P101SR00","plasma_P101SR00","plasma_P105SR00","plasma_P106SR00"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#FF450080","#FF450080","#FF450080","#FF450080","#FF450080","#FF450080","#00CD6680","#00CD6680","#00CD6680","#00CD6680","#00CD6680","#00CD6680","#00B2EE80","#00B2EE80","#00B2EE80","#00B2EE80","#00B2EE80","#00B2EE80"),do.return = T ,pt.size = 2,do.label = F,label.size=0,plot.title="merge_mutation")

pdf(file = paste("106.pdf"),width=14,height = 8)
p1<-TSNEPlot(object = SampleAlignment106,group.by = "cellident",do.return = T ,pt.size = 1,do.label = T,label.size=0,plot.title="P106",no.legend =F)
p2<-TSNEPlot(object = SampleAlignment,group.by = "origident",do.return = T ,pt.size = 1,do.label = T,label.size=3,cells.use=celluse105,plot.title="P105_tumor",no.legend =F)
p2<-TSNEPlot(object = SampleAlignment,group.by = "origident",plot.order=c("stemCell_1","stemCell_2","stemCell_3","pro_1","CD14Mono","cDC","B_1","B_2","CD16Mono_1","CD16Mono_2","CD4 memo T","CD4 Treg","CD4T","CD8 memo T","CD8T","Eryth","Granu_1","Granu_2","MK","N", "NK","NKT","pDC","pro_2","pro_3","stemCell_4"),colors.use=c("#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#00B2EE","#B452CD","#43CD80","#1874CD","#FFC125","#FF4040"),cells.use=celluse105,do.return = T ,pt.size = 1,do.label = TRUE,label.size=0,plot.title="P105_tumor",no.legend =F)
plot_grid(p1,p2)
dev.off()


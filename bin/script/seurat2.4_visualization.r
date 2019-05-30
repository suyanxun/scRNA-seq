
#===============================================
#计算突变在的群体
SampleAlignmentSeurat<-readRDS(file = paste(patientID,"_",dim,"SampleAlignmentSeurat.rds",sep=""))

mutationNo<-table(SampleAlignmentSeurat@meta.data$orig.ident,SampleAlignmentSeurat@meta.data$finalCellMutaionStatus_PTPN11)
onlyMut<-subset(SampleAlignmentSeurat@meta.data,finalCellMutaionStatus_PTPN11=="Mut")
#计算每个群体在每个样本中的突变数,用于确定肿瘤群体(D0>2)
mutationNoInClusterStim<-table(onlyMut$orig.ident,onlyMut$stim)

#计算比例
clusterNo<-table(SampleAlignmentSeurat@ident,SampleAlignmentSeurat@meta.data$stim)
clusterRatio<-prop.table(clusterNo,2)*100

#计算比例（只看肿瘤群体）
#blast<-c("stemCell_1","pro_1","stemCell_2","cDC","pDC","stemCell_3","pro_2","pro_3")  #FLT3
blast<-c("stemCell_1","CD14Mono","stemCell_2","cDC","stemCell_3","pro_1") #D0突变，与BM涂片结果大致对应
colSums(clusterRatio[blast,])
#P105SR00-BM P105SR00-PB P105SR01-PB P105SR04-PB P105SR26-BM P105SR26-PB 
#  41.915246   26.894437   25.109451    7.649214   15.339117    4.349186

#合并,输出为表格
library(tidyr)
mutationNo<-as.data.frame(mutationNo)
mutationNoInClusterStim<-as.data.frame(mutationNoInClusterStim)
mutationNo<-spread(mutationNo,Var2,Freq)
colnames(mutationNo)[1]<-"Cluster"
mutationNoInClusterStim<-spread(mutationNoInClusterStim,Var2,Freq)
colnames(mutationNoInClusterStim)<-c("Cluster",paste(colnames(mutationNoInClusterStim)[2:7],"MutNo",sep="_"))

clusterNo<-as.data.frame(clusterNo)
clusterNo<-spread(clusterNo,Var2,Freq)
colnames(clusterNo)<-c("Cluster",paste(colnames(clusterNo)[2:7],"clusterNo",sep="_"))
rownames(clusterNo)<-clusterNo[,1]
clusterRatio<-as.data.frame(clusterRatio)
clusterRatio<-spread(clusterRatio,Var2,Freq)
colnames(clusterRatio)<-c("Cluster",paste(colnames(clusterRatio)[2:7],"clusterRatio",sep="_"))
rownames(clusterRatio)<-clusterRatio[,1]

count<-merge(mutationNoInClusterStim,mutationNo,by="Cluster") 
count<-merge(count,clusterNo,by="Cluster")
count<-merge(count,clusterRatio,by="Cluster")
write.table(count,paste(patientID,"_",dim,"Samplemerge_Count_",res,".txt",sep=""),sep="\t",quote = FALSE,row.names=FALSE)

#=================================ratio changes of each cluster=========================
#show the ratio changes of each cluster 
Freq<-table(SampleAlignmentSeurat@meta.data$stim,SampleAlignmentSeurat@meta.data$orig.ident,SampleAlignmentSeurat@meta.data$Phase)
#a<-table(SampleAlignmentSeurat@meta.data$stim,SampleAlignmentSeurat@meta.data$Phase)
#addmargins(a)
Freq<-as.data.frame(Freq)
colnames(Freq)<-c("SampleID","Cluster","Phase","Count")
Cluster_Percet = ddply(Freq,"SampleID",transform,percent = Count/sum(Count) * 100)
Cluster_Percet$SampleID<-factor(Cluster_Percet$SampleID,levels=sampleName)
Cluster_Percet$Phase<-factor(Cluster_Percet$Phase,levels=c("G2M","S","G1"))

colourCount = length(unique(Cluster_Percet$Cluster))
getPalette = colorRampPalette(colors = brewer.pal(9, "Set1"))

pdf(file = paste(patientID,"_",dim,"Samplemerge_ClusterPercetAndCellCycleChanges_res",res,".pdf",sep=""),width=45,height = 25)
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
TumorclusterRatioInAllCluster<-Cluster_Percet[which(Cluster_Percet$Cluster %in% blast),]
pdf(file = paste(patientID,"_",dim,"Samplemerge_ClusterPercetAndCellCycleChanges_TumorLike_res",res,".pdf",sep=""),width=45,height = 25)
p1<-ggplot(data=TumorclusterRatioInAllCluster, aes(x=SampleID, y=percent, fill=Cluster))+
  geom_bar(stat="identity")+
  #	geom_text(aes(label = Count), vjust = 1.5, colour = "black", position = position_dodge(.9), size = 2)+
  labs(title=paste(dim,"clusterPercentChanges_res_",res,sep=""))+
  theme(axis.text.x = element_text(angle=90,size = 30),axis.text.y = element_text(size = 40),text = element_text(size = 40),legend.position = "top")+
#  scale_fill_brewer(palette = "Set3",direction = -1)+
  scale_fill_manual(values = getPalette(colourCount)) 
#  facet_grid(.~Cluster)
plot(p1)
dev.off()

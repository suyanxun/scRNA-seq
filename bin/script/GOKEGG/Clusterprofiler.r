setwd("F://scRNAseq/Result/IntegratingMultiSample/SixSampleAlignment_20190107/")

#source("http://bioconductor.org/biocLite.R")
#biocLite('clusterProfiler')
#biocLite("org.Hs.eg.db")

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

df<-read.table(file="./geneExpressInRemovedCell.txt",header=F,sep="\t")
colnames(df)<-"gene"
gene<-as.character(df$gene)


#table(is.na(mygene))

for (i in 1:(length(gene))) {
  temp<-strsplit(gene[i],"\\.")
  gene[i]<-temp[[1]][1]
}

uniqueGene<-unique(gene)
mygene<-select(org.Hs.eg.db,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL",keys=uniqueGene)
table(is.na(mygene)) #NA number
ego<-enrichGO(OrgDb = "org.Hs.eg.db",gene=mygene$ENTREZID,ont="BP",pvalueCutoff=0.01,readable=TRUE)

write.table(ego,"GO-enrich_cluster.tsv",sep="\t",row.names = F)
pathway<-as.data.frame(ego)
pathway = pathway[1:30,]
pathway = pathway[order(pathway$pvalue,decreasing=T),]
tt <- factor(pathway$Description, levels=unique(pathway$Description))
pp = ggplot(pathway,aes(-1*log10(p.adjust),tt))
pbubble = pp +geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_colour_gradient(low="blue",high="red") +
  theme_bw() +scale_size(range=c(2, 10)) +
  labs(title = "GO-enrich_cluster")+
  ylab("")+xlim(min(-1*log10(pathway$p.adjust))*0.9, max(-1*log10(pathway$p.adjust))*1.1)+
  theme(axis.text.y = element_text(size = 25, family = "Helvetica", color = "black",  angle = 0),axis.text.x = element_text(size = 15, family = "Helvetica", color = "black", angle = 0))+
  theme(legend.title = element_text(color="black", size=20, family = "Helvetica"))+
  theme(legend.text = element_text(color="azure4", size = 15,  family = "Helvetica"))+
  theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
  theme(legend.key.size=unit(1.1,'cm'))
ggsave(file="GO-enrich_cluster.pdf", plot=pbubble, width=20, height=15)

# KEGG enrichment analysis  
ekk<-enrichKEGG(gene=mygene$ENTREZID,organism = "hsa",pvalueCutoff=0.05)
write.table("KEGG-enrich_cluster.tsv",row.names = F)
pathway<-as.data.frame(ekk)
pathway = pathway[1:30,]
pathway = pathway[order(pathway$pvalue,decreasing=T),]
tt <- factor(pathway$Description, levels=unique(pathway$Description))
pp = ggplot(pathway,aes(-1*log10(p.adjust),tt))
pbubble = pp +geom_point(aes(size=Count,color=-1*log10(p.adjust)))+
  scale_colour_gradient(low="blue",high="red") +theme_bw() +
  scale_size(range=c(4, 12)) +
  labs(title = "KEGG-enrich")+ylab("")+
  theme(axis.text.y = element_text(size = 25, family = "Helvetica", color = "black", angle = 0),axis.text.x = element_text(size = 15, family = "Helvetica", color = "black",  angle = 0))+
  theme(legend.title = element_text(color="black", size=20, family = "Helvetica"))+
  theme(legend.text = element_text(color="azure4", size = 15,  family = "Helvetica"))+
  theme(axis.title = element_text(color="black", size=16, family = "Helvetica"))+
  theme(legend.key.size=unit(1.1,'cm'))
ggsave(file="KEGG-enrich_cluster.pdf", plot=pbubble, width=15, height=15)






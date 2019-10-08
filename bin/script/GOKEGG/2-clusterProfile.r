#根据gene_list进行GO及KEGG富集分析

#/pnas/wangqf_group/yuxx/software/miniconda3/bin/R

library(ini, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))
library(optparse, lib.loc=c("/pnas/wangqf_group/suyx/PMO/bin/R3.5_package"))

usage <- function(){ #使用说明
	print_help(opt_parser)
}

option_list = list(
	make_option( c("-g", "--gene_list"), type="character", default=NULL, 
				help="[Required] geneList文件，一列gene symbol", metavar = "character"),
	make_option(c("-o", "--outdir"), type="character", default=NULL, 
				help="[Required] 结果输出目录", metavar="character"),
	make_option(c("-s", "--sample"), type="character", default=NULL, 
				help="[Optional] 样本名", metavar="character")
);
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$g) || is.null(opt$o)){
	usage()
	q()
}
if(!file.exists( opt$o )) {
	dir.create(opt$o)
}
if (is.null(opt$s)){
	opt$s = basename( opt$g )
}

library(DOSE)
library(org.Hs.eg.db)
#library(topGO)
library(clusterProfiler)
#library(pathview)
library(ggplot2)

gene <- read.table( opt$g, header=FALSE)
gene$V1 <- as.character(gene$V1)
gene_trans = bitr(gene$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
data(geneList, package="DOSE")
ego_ALL <- enrichGO(gene = gene_trans$ENTREZID,universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = TRUE)
write.table(as.data.frame(ego_ALL),paste(opt$o,"/",opt$s,"_ALL-enrichGo.csv",sep=""),row.names =FALSE, sep = "\t", quote = FALSE)
plots = dotplot(ego_ALL,title= paste(opt$s," EnrichmentGO",sep = ""))
ggsave( paste(opt$o,"/",opt$s,"_ALL-EnrichmentGO_dot.png",sep = ""), plots, height = 10, width = 10 )
plots = barplot(ego_ALL, showCategory=20,title="EnrichmentGO_ALL")
ggsave( paste(opt$o,"/",opt$s,"_ALL-EnrichmentGO_bar.png",sep = ""), plots, height = 10, width = 10 )

ekk = enrichKEGG(gene = gene_trans$ENTREZID, organism = "hsa", universe = names(geneList), pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, use_internal_data = FALSE)
write.table(as.data.frame(ekk@result), file=paste(opt$o,"/",opt$s,"_ALL-enrichKegg.csv",sep=""), sep = "\t", quote = FALSE)


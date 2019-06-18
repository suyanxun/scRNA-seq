library(dplyr)
library(Seurat, lib.loc=c("/asnas/wangqf_group/jiangsht/Software/R"))
library(ggplot2)
library(cowplot)


#dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/1_CellRanger_rawdata/P101SR00-PB/cellrangerRes/P101SR00-PB/outs/filtered_feature_bc_matrix"
#seurat_object = as.matrix(Read10X(data.dir = dir))

dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/all_Dim20_Res2.4_mark_mutation_and_celltype_FAB.rds"
seurat_object = readRDS(dir)
gene = row(seurat_object[["integrated"]]@data)
temp = GetAssayData( object = seurat_object, slot = "counts")

express_matrix = as.matrix(temp[gene,])
seurat_object = express_matrix
result_row = dim(seurat_object)[1] + 4
result_col = dim(seurat_object)[2] + 2


result = matrix( data=rep("",times=result_row*result_col),nrow=result_row,ncol=result_col)
result[5:result_row,3:result_col] = seurat_object
result[1,1:7] = c("CEF","1","1","1",dim(seurat_object)[1], dim(seurat_object)[2], 0)
result[2,1:2] = c("header","AML")
result[3,] = c("","CellID",colnames(seurat_object))
result[4,1] = "gene"
result[5:result_row,1] = rownames(seurat_object)
write.table(result,"/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/all/backspin/express_matrix.cef",sep="\t",quote = FALSE,row.names=FALSE,col.names=FALSE)

###使用python运行backSPIN分析####


###绘制热图查看###
library( pheatmap )
library( ggplot2 )
result = as.matrix(read.table( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P101/backspin/backspin_result.cef", sep = "\t" ))
express_matrix = result[ (dim(result)[1] - as.numeric(result[1,5])+1 ):dim(result)[1], (dim(result)[2] - as.numeric(result[1,6])+1):dim(result)[2] ]
express_matrix = apply( express_matrix, 2, as.numeric )
colnames(express_matrix) = result[ as.numeric(result[1,2])+2, (as.numeric(result[1,3])+2):dim(result)[2]]
rownames(express_matrix) = result[ (as.numeric(result[1,2]) + as.numeric(result[1,4]) + 3):dim(result)[1],1]
cluster = result[(as.numeric(result[1,2]) + as.numeric(result[1,4]) + 1),(as.numeric(result[1,3])+2):dim(result)[2]]
names(cluster) = colnames(express_matrix)
express_matrix = log( express_matrix+1 )
annotation_col = data.frame(Cluster = factor(cluster) )
rownames(annotation_col) = colnames(express_matrix)
#plots = pheatmap(express_matrix, cluster_rows = FALSE, cluster_cols = FALSE,show_rownames = FALSE, show_colnames = FALSE,annotation_col = annotation_col )
plots = pheatmap(express_matrix, cluster_cols = FALSE,show_rownames = FALSE, show_colnames = FALSE,annotation_col = annotation_col )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P101/backspin/a.png", plots )


###增加突变信息
mutation = c( "CTAAGACGTCTTCGTC","GAATGAAAGCACGCCT", "AAGCCGCTCCGCGGTA", "ACATGGTCAACACGCC", "ATTATCCGTTGTCGCG", "CAACCAAAGAAGATTC", "CGGACGTTCCGTTGCT", "CTGAAGTGTTACGACT", "GACCAATGTCTAAACC", "GACGTTAAGACGCAAC", "GATCGCGCACATCTTT", "GGACATTTCTTAACCT", "GGGACCTCACGCGAAA", "GTGCAGCTCAGCATGT", "TACTCATCACGAGGTA" )
mutation_type = table(cluster[mutation])

mutation_number = rep("0",length( cluster ))
for (i in names(mutation_type)){
	mutation_number[cluster == i] = as.character( mutation_type[i] )
}

annotation_col = data.frame(Cluster = factor(cluster), mutation_number =  factor(mutation_number) )
rownames(annotation_col) = colnames(express_matrix)
#plots = pheatmap(express_matrix, cluster_rows = FALSE, cluster_cols = FALSE,show_rownames = FALSE, show_colnames = FALSE,annotation_col = annotation_col )
plots = pheatmap(express_matrix, cluster_cols = FALSE,show_rownames = FALSE, show_colnames = FALSE,annotation_col = annotation_col, scale="row" )
ggsave( "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/P101/backspin/b.png", plots )



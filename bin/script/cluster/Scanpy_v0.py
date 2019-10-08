import scanpy as sc
import numpy as np
import pandas as pd

h5_file = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/1_CellRanger_rawdata/P110SR00-BM/cellrangerRes/P110SR00-BM/outs/filtered_feature_bc_matrix.h5"

file_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/1_CellRanger_rawdata/P110SR00-BM/cellrangerRes/P110SR00-BM/outs/filtered_feature_bc_matrix"

adata = sc.read_10x_mtx( file_dir )
adata = sc.read_10x_h5( h5_file )


adata.var_names_make_unique()
sc.pl.highest_expr_genes(adata, n_top=20) #可视化表达占比最高的20个基因


#基础过滤：
sc.pp.filter_cells(adata, min_genes=200)   # 去除表达基因200以下的细胞
sc.pp.filter_genes(adata, min_cells=3)     # 去除在3个细胞以下表达的基因

#质控
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

#小提琴图 及 散点图
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

#过滤基因及线粒体百分比
adata = adata[adata.obs['n_genes'] < 4000, :]
adata = adata[adata.obs['percent_mito'] < 0.3, :]

#数据标准化
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata

#识别差异基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

#保留差异分析基因用于后续分析
adata = adata[:, adata.var['highly_variable']]

#线性回归
sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
sc.pp.scale(adata, max_value=10)

#主成分分析
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write("/pnas/wangqf_group/suyx/PMO/test/scanpy/pca_results.h5ad")

#聚类分析
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)


sc.pl.umap(adata, color=['louvain'])
adata.write("/pnas/wangqf_group/suyx/PMO/test/scanpy/umap.h5ad")



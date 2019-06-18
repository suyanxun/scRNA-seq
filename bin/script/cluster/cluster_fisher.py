'''
根据fisher检验进行分类
'''
import os
import sys
import argparse
import logging
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency

#所有的理论数T≥5并且总样本量n≥40,用Pearson卡方进行检验.
#如果理论数T＜5但T≥1,并且n≥40,用连续性校正的卡方进行检验.
#如果有理论数T＜1或n＜40,则用Fisher’s检验.


__author__ = 'suyanxun'
__mail__ = 'suyx1006@163.com'

def load_healthy_marker( marker_file ):
	'''
	读取正常人样本分类的singernature gene
	'''
	reasult = {}
	with open( marker_file, "r" ) as infile:
		for line in infile:
			gene, cluster = line.strip().split("\t")
			if cluster in reasult:
				reasult[cluster].append( gene )
			else:
				reasult[cluster] = [gene]
	return( reasult )

def load_seurat_marker( marker_file ):
	'''
	读取seurat分类结果
	'''
	result = {}
	n = 0
	current_cluster = ""
	all_gene = []
	with open( marker_file, 'r' ) as infile:
		for line in infile:
			if line.startswith( "p_val" ):continue
			*temp, cluster, gene = line.strip().split( "\t" )
			if cluster != current_cluster:
				n = 0
			n += 1
			if n <= 50:
				if cluster not in result: result[ cluster ] = [gene]
				else: result[ cluster ].append( gene )
				all_gene.append( gene )
			else:
				all_gene.append( gene )
	return result	

def test_cluster( healthy_marker, seurat_cluster_gene ):
	result = "unknown"
	p = 1
	total_gene = 1500
	for cluster in healthy_marker:
		healthy_gene = healthy_marker[ cluster ]
		all_in = len( [ i for i in healthy_gene if i in seurat_cluster_gene ])
		healthy_only = len( list(set(healthy_gene).difference(set(seurat_cluster_gene))) )
		seurat_only = len( list(set(seurat_cluster_gene).difference(set(healthy_gene))) )
		if all_in < 1:
			_, p_val = fisher_exact( [[ all_in, healthy_only],[seurat_only,total_gene]] )
		elif all_in >= 5 and healthy_only >= 5 and seurat_only >= 5:
			_, p_val, *temp = chi2_contingency([[ all_in, healthy_only],[seurat_only,total_gene]],False)
		else:
			_, p_val, *temp = chi2_contingency([[ all_in, healthy_only],[seurat_only,total_gene]],True)
		if p_val < p and p_val< 0.01:
			result = cluster
			p = p_val
	return result

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--infile',help='seurat marker result file',dest='infile',required=True)
	parser.add_argument('-o','--outfile',help='output file',dest='outfile',required=True)
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	logging.info("开始分析")
	marker_file = "/pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/config/CellSignatureTumorDerived.txt"
	healthy_marker = load_healthy_marker( marker_file )
	seurat_cluster = load_seurat_marker( args.infile )
	result = {}
	for cluster in seurat_cluster:
		seurat_cluster_gene = seurat_cluster[cluster]
		seurat_cell_type = test_cluster( healthy_marker, seurat_cluster_gene )
		if seurat_cell_type in result:
			result[ seurat_cell_type ].append( cluster )
		else:
			result[ seurat_cell_type ] = [ cluster ]
	with open( args.outfile, "w" ) as outfile:
		for seurat_cell_type in result:
			outfile.write( "{0}={1}\n".format( seurat_cell_type, ",".join( result[seurat_cell_type] ) ) )

if __name__=="__main__":
	main()

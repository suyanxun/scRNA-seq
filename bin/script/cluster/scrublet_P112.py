# # -*- coding:utf-8 -*-
# '''
# 使用scrublet模块鉴定双细胞
# '''
import os
import sys
import argparse
import skimage.filters
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import logging
sys.path.append("/pnas/wangqf_group/suyx/PMO/bin/python3_package")
#sys.path.append("/asnas/wangqf_group/jiangsht/Software/anaconda3/bin")
import scrublet as scr
import scipy.io


__author__ = 'suyanxun'
__mail__ = 'suyx1006@163.com'

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-m','--mtx',help='cellranger分析结果中的matrix.mtx',dest='mtx',required=True)
	parser.add_argument('-f','--feature',help='cellranger分析结果中的feature.csv或genes.csv',dest='feature',required=True)
	parser.add_argument('-o','--outdir',help='结果输出目录',dest='outdir',required=True)
	parser.add_argument('-s','--sampleName',help='样本名',dest='sampleName',required=True)
	parser.add_argument('-e','--expectedDoubletRate',help='细胞结团率',dest='expectedDoubletRate',type=float, required=True)
	parser.add_argument('-p','--pc',help='PC值',dest='pc',type = int, default = 30)
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	logging.info("开始分析")

	# plt.rcParams['font.family'] = 'sans-serif'
	# plt.rcParams['font.sans-serif'] = 'Arial'
	# plt.rc('font', size=14)
	# plt.rcParams['pdf.fonttype'] = 42

	pc = args.pc
	expectedDoubletRate = args.expectedDoubletRate
	sampleName= args.sampleName

	#Load counts matrix and gene list
	counts_matrix = scipy.io.mmread( args.mtx ).T.tocsc()
	genes = np.array(scr.load_genes( args.feature, delimiter='\t', column=1))
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expectedDoubletRate)
	
	#Run the default pipeline
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=pc)
	#Plot doublet score histograms for observed transcriptomes and simulated doublets
	scrub.call_doublets(threshold=0.1)
	scrub.plot_histogram()
	plt.savefig( args.outdir + "/" + sampleName + '_Scrublet_Histogram.png')
	scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
	scrub.plot_embedding('UMAP', order_points=True)
	plt.savefig( args.outdir + "/" + sampleName + '_Scrublet_UMAP.png')

	#output the log file
	#expected doublet rate、 detected doublet rate、 doublet threshold、overall doublet rate
	logging.info( "patientID\texpected_doublet_rate\tdetected_doublet_rate\toverall_doublet_rate\tthreshold\tPC\n" )
	logging.info("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n" %(sampleName,scrub.expected_doublet_rate,scrub.detected_doublet_rate_,scrub.overall_doublet_rate_,scrub.threshold_,pc))

	#output the doublet status of every single cell
	with open(args.outdir + "/" + sampleName + ".predictDoublet_scrublet.txt", "w+") as fo:
		for i in scrub.predicted_doublets_ :
			fo.write("%s\n" %(i))

if __name__=="__main__":
	main()


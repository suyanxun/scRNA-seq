'''
生成spring软件目录
'''
import os
import sys
import argparse
import logging
import numpy as np

sys.path.append("/pnas/wangqf_group/suyx/PMO/softwore/SPRING")
from preprocessing_python import *


__author__ = 'suyanxun'
__mail__ = 'suyx1006@163.com'

def load_express_matrix( express_matrix_file ):
	'''
	表达矩阵，每行一个细胞，每列一个gene，
	'''
	express_matrix = []
	gene_list = []
	cell_ID = []
	with open( express_matrix_file, "r" ) as infile:
		line_num = 0
		for line in infile:
			line_num += 1
			if line_num == 1:
				gene_list = line.strip().split( "\t" )
				continue
			tmp = line.strip().split( "\t" )
			cell_ID.append( tmp[0] )
			express_matrix.append( [float(i) for i in tmp[1:] ] )
	express_matrix = np.array(express_matrix)
	return( express_matrix, cell_ID, gene_list )

def load_attributes( attributes_file ):
	'''
	每行为一个属性，第一列为属性名。每个细胞一列
	'''
	attributes = {}
	with open( attributes_file, 'r' ) as infile:
		for line in infile:
			tmp = line.strip().split( "\t" )
			attributes[ tmp[0] ] = tmp[1:]
	return( attributes )

def generate_project( express_matrix, gene_list, attributes, cell_ID, project_dir ):
	
	custom_colors = { "uniform":[1]*len(cell_ID) }
	cell_groupings = attributes
	Epca = get_PCA(Zscore(express_matrix),20)
	
	D = get_distance_matrix( Epca )
	E = express_matrix
	save_spring_dir(E,D,5,gene_list,project_dir, cell_groupings=cell_groupings, custom_colors=custom_colors, coarse_grain_X=1)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-e', '--express_matrix', help='express matrix file', dest='express_matrix', required=True)
	parser.add_argument('-a', '--attributes', help='attributes file', dest='attributes', required=True)
	parser.add_argument('-o', '--outdir', help='project dir', dest='outdir', required=True)
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	
	softwore_dir = "/pnas/wangqf_group/suyx/PMO/softwore/SPRING"
	
	express_matrix, cell_ID, gene_list = load_express_matrix( args.express_matrix )
	attributes = load_attributes( args.attributes )
	generate_project( express_matrix, gene_list, attributes, cell_ID, args.outdir + "/outs" )
	os.system( "cp -r {0}/springViewer.html {1} ".format( softwore_dir, args.outdir ) )
	os.system( "cp -r {0}/scripts {1} ".format( softwore_dir, args.outdir ) )
	os.system( "cp -r {0}/stuff {1} ".format( softwore_dir, args.outdir ) )
	

if __name__=="__main__":
	main()

'''
使用phenogra包进行分析
'''
import os
import sys
import argparse
import logging
import phenograph
import numpy

__author__ = 'suyanxun'
__mail__ = 'suyx1006@163.com'

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-s','--seq',help='fasta file',dest='seq',required=True)
	parser.add_argument('-s','--seq',help='fasta file',dest='seq',default = None)
	parser.add_argument('-s','--seq',help='fasta file',action='store_true')
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	
	file = open("/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/PublicData/cell/GSE116256/suppl/GSM3587923_AML1012-D0.dem.txt")
	data = file.readlines()
	file.close()
	temp = data[:]
	cellID = data[0].strip().split("\t")[1:]
	gene = [ i.strip().split("\t")[0] for i in data[1:] ]
	counts = [ i.strip().split("\t")[1:] for i in data[1:] ]
	counts = [ [ int(j) for j in i ] for i in counts ]
	counts_np = numpy.array( counts )
	counts_np = counts_np.T
	communities, graph, Q = phenograph.cluster(counts_np)

if __name__=="__main__":
	main()

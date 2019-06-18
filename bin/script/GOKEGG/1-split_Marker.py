'''
分割marker结果，输出各个cluster中对应的满足过滤条件的gene list
'''
import os
import sys
import argparse
import logging

__author__ = 'suyanxun'
__mail__ = 'suyx1006@163.com'

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:	{0}	mail:	{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--infile',help='marker结果文件',dest='infile',required=True)
	parser.add_argument('-p','--p_adj',help='矫正后p值筛选阈值',dest='p', type = float, default = 0.01)
	parser.add_argument('-f','--fc',help='foldChange值（取log之后）筛选阈值(绝对值)',dest='fc',type = float,default = 0)
	parser.add_argument('-o','--outdir',help='结果输出目录',dest='outdir',required=True)
	
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	result = {}
	with open( args.infile, "r" ) as infile:
		for line in infile:
			if "avg_logFC" in line:continue
			tmp = line.strip().split( "\t" )
			cluster = tmp[-2]
			gene = tmp[-1]
			p_val_adj = float( tmp[-3] )
			avg_logFC = float( tmp[2] )
			if p_val_adj < args.p and abs(avg_logFC) > args.fc:
				type = "negative" if avg_logFC < 0 else "positive"
				if cluster in result:
					if type in result[cluster]:
						result[cluster][type].append( gene )
					else:
						result[cluster][type] = [gene]
				else:
					result[cluster] = {type:[gene]}
	for cluster in result:
		all_gene = []
		for type in ["negative", "positive"]:
			if type in result[cluster]:
				with open( "{0}/{1}_{2}.txt".format(args.outdir,cluster,type) , "w") as outfile:
					outfile.write( "{0}\n".format( "\n".join(result[cluster][type]) ) )
					all_gene.extend( result[cluster][type] )
		with open( "{0}/{1}_{2}.txt".format(args.outdir,cluster,"all") , "w") as outfile:
			outfile.write( "{0}\n".format( "\n".join(all_gene) ) )

if __name__=="__main__":
	main()

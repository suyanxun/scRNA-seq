'''
筛选Annovar注释之后的结果，筛选标准：
1. 位于正常染色体上
2. 覆盖度在50以上
3. 位于外显子区域上（需要参数控制）
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
	parser.add_argument( '-i', '--input', help='输入文件', dest='input', required=True )
	parser.add_argument( '-o', '--output', help='输出文件', dest='output', required=True )
	parser.add_argument( '-d', '--depth', help='最低测序深度', dest='depth', type=int, default = 50 )
	parser.add_argument( '-e', '--exon', help='是否只保留外显子区域位点', action='store_true' )
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	
	#stat = [ "sample", "all_siteNum", "E0>50", "EC>50", "exonic", "filter" ]
	stat = [ os.path.basename( args.input ), 0, 0, 0, 0, 0 ]
	tmp = [ i+1 for i in range(22)]
	tmp.extend( ["X","Y"] )
	chrs = [ "chr"+str(i) for i in tmp ]
	with open( args.input, "r") as input, open( args.output, "w" ) as output:
		for line in input:
			if line.startswith("CHROM"):
				output.write(line)
			else:
				stat[1] += 1
				line_info = line.strip().split("\t")
				chr = line_info[0]; func = line_info[5]
				E0 = line_info[11]; EC = line_info[14]
				if chr in chrs and int(E0) > args.depth and int(EC) > args.depth:
					if args.exon:
						if "exonic" in func:
							output.write( line )
							stat[5] += 1
					else:
						output.write( line )
						stat[5] += 1
				if int(E0) <= args.depth:stat[2] += 1
				if int(EC) <= args.depth:stat[3] += 1
				if "exonic" in func:stat[4] += 1
	with open( args.output + ".stat.txt", "w" ) as output:
		output.write( "\t".join(["sample", "all_siteNum", "E0<50", "EC<50", "exonic", "filter"]) + "\n")
		output.write( "\t".join( [str(i) for i in stat] ))
	

if __name__=="__main__":
	main()

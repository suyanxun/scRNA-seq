'''
统计cb_sniffer输出结果: counts_CB.tsv文件
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
	parser.add_argument( '-i', '--infiles', help='输入文件', dest='infiles', required=True, nargs="+" )
	parser.add_argument( '-o', '--outfile', help='输出文件', dest='outfile', required=True )
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	
	result = {}
	output = open( args.outfile, "w" )
	output.write( "Sample\tWT_cell\tMU_cell\tMix_cell\tGene_num\tCell_num\n" )
	for infile in args.infiles:
		sample_name = os.path.basename( infile ).replace( "_counts_CB.tsv", "" )
		output.write( sample_name + "\t" )
		if sample_name not in result:
			result[sample_name] = [0,0,0]
		cell = []
		gene = []
		with open(infile, 'r') as input:
			for line in input:
				if line.startswith("chrm"):continue
				*_, gene_name,cell_ID,ref_count,alt_count,_ = line.strip().split("\t")
				ref_count = int(ref_count); alt_count = int(alt_count)
				gene.append(gene_name); cell.append(cell_ID)
				if ref_count > alt_count and alt_count == 0:
					result[sample_name][0] += 1
				elif alt_count > ref_count and ref_count == 0:
					result[sample_name][1] += 1
				else:
					result[sample_name][2] += 1
		gene_num = len(set(gene))
		cell_num = len(set(cell))
		result[sample_name].extend( [gene_num, cell_num] )
		result[sample_name] = [ str(i) for i in result[sample_name] ]
		output.write( "\t".join( result[sample_name] ) + "\n")
	output.close()
 
 
if __name__=="__main__":
	main()

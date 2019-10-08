'''
根据突变筛选结果及单细胞突变筛选结果生成点亮突变细胞的结果文件
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
	parser.add_argument('-i1','--infile1',help='靶向数据突变筛选结果',dest='infile1',required=True)
	parser.add_argument('-i2','--infile2',help='单细胞突变筛选结果',dest='infile2',required=True)
	parser.add_argument('-o','--outfile',help='结果输出目录',dest='outfile',required=True)
	parser.add_argument('-s','--sample',help='样本名',dest='sample',required=True)
	args=parser.parse_args()

	logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s - %(message)s")
	
	mutation_site = []
	with open( args.infile1, "r" ) as infile1:
		for line in infile1:
			if line.startswith("CHROM"):continue
			tmp = "_".join(line.strip().split("\t")[0:3]).replace( "chr", "" )
			mutation_site.append(tmp)
	with open( args.infile2, "r" ) as infile2, open( args.outfile, "w" ) as outfile:
		exist_site = []
		outfile.write( "\t".join([ "allTransNo", "patientIDCB", "patientID", "gene", "position",  "transMutNo","transWTNo", "transFurther2confirmNo", "cellMutaionStatus", "finalCellMutaionStatus", "coverInMutPos", "allTransNoFilter" ]) + "\n" )
		for line in infile2:
			if line.startswith("chrm"):continue
			tmp = line.strip().split("\t")
			if "_".join(tmp[0:3]) in mutation_site:
				patientIDCB = "_".join( [ args.sample, tmp[7].split("-")[0] ] )
				if patientIDCB in exist_site:continue
				else: exist_site.append(patientIDCB)
				patientID = args.sample
				gene = tmp[6]
				position = tmp[1]
				allTransNo = tmp[10]
				transMutNo = tmp[9]
				transWTNo =  tmp[8]
				transFurther2confirmNo = "0"
				cellMutaionStatus = "Mut" if int(transMutNo) > int(transWTNo) else "WT"
				finalCellMutaionStatus = cellMutaionStatus
				coverInMutPos = "-"
				allTransNoFilter = allTransNo
				outfile.write( "\t".join([ patientIDCB, allTransNo, patientIDCB, patientID, gene, position, transMutNo, transWTNo, transFurther2confirmNo, cellMutaionStatus, finalCellMutaionStatus, coverInMutPos, allTransNoFilter ]) + "\n" )
	
	

	 
if __name__=="__main__":
	main()

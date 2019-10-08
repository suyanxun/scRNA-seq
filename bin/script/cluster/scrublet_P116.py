# -*- coding:utf-8 -*-

#/software/biosoft/software/python/python3/bin/python
import sys
import scrublet as scr
import scipy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

pc=30
#expectedDoubletRate=0.0264
#sampleName="P105SR00-BM"

#expectedDoubletRate=0.0138
#sampleName="P106SR00-BM"

#expectedDoubletRate=0.0167
#sampleName="P107SR00-BM"

#expectedDoubletRate=0.0101
#sampleName="P108SR00-BM"


#expectedDoubletRate=0.0147
#sampleName="P110SR00-BM"

#expectedDoubletRate=0.0239
#sampleName="P111SR00-BM"

#expectedDoubletRate=0.0192
#sampleName="P114SR00-BM"

#expectedDoubletRate=0.0077
#sampleName="P115SR00-BM"

#expectedDoubletRate=0.0258
#sampleName="P116SR00-BM"

#expectedDoubletRate=0.0342
#sampleName="P117SR00-BM"

expectedDoubletRate=0.0092
sampleName="P119SR00-BM"

wkdir="/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/scrublet_result/"
#用于改变当前工作目录到指定的路径
os.chdir(wkdir)

input_dir = "/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/2_seurat/diagnosis/" + sampleName + "/"
#Load counts matrix and gene list
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
genes = np.array(scr.load_genes(input_dir + '/features.tsv', delimiter='\t', column=1))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expectedDoubletRate)
#Run the default pipeline
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=pc)
#Plot doublet score histograms for observed transcriptomes and simulated doublets
#scrub.call_doublets(threshold=0.3)   #"P105SR00-BM"
#scrub.call_doublets(threshold=0.3)   #"P106SR00-BM"
#scrub.call_doublets(threshold=0.2)   #"P107SR00-BM"
#scrub.call_doublets(threshold=0.1)   #"P108SR00-BM"
#scrub.call_doublets(threshold=0.5)   #"P110SR00-BM"
#scrub.call_doublets(threshold=0.4)   #"P111SR00-BM"
#scrub.call_doublets(threshold=0.32)   #"P114SR00-BM"
#scrub.call_doublets(threshold=0.21)   #"P115SR00-BM"
#scrub.call_doublets(threshold=0.23)   #"P116SR00-BM"
#scrub.call_doublets(threshold=0.5)   #"P117SR00-BM"
scrub.call_doublets(threshold=0.14)   #"P119SR00-BM"


scrub.plot_histogram()
plt.savefig('./'+sampleName+'_Scrublet_Histogram.png')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.savefig('./'+sampleName+'_Scrublet_UMAP.png')

#output the log file
#expected doublet rate、 detected doublet rate、 doublet threshold、overall doublet rate
if not os.path.isfile(wkdir+"/"+sampleName+".predictDoublet_scrublet.log"):
    log = open(wkdir+"/"+sampleName+".predictDoublet_scrublet.log", "a+")
    log.write("patientID\texpected_doublet_rate\tdetected_doublet_rate\toverall_doublet_rate\tthreshold\tPC\n")
    log.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n" %(sampleName,scrub.expected_doublet_rate,scrub.detected_doublet_rate_,scrub.overall_doublet_rate_,scrub.threshold_,pc))
else:
    log = open(wkdir+"/"+sampleName+".predictDoublet_scrublet.log", "a+")
    log.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%s\n" %(sampleName,scrub.expected_doublet_rate,scrub.detected_doublet_rate_,scrub.overall_doublet_rate_,scrub.threshold_,pc))

log.close()

#output the doublet status of every single cell
fo = open("./"+sampleName+".predictDoublet_scrublet.txt", "w+")
for i in scrub.predicted_doublets_ :
    fo.write("%s\n" %(i))

fo.close()

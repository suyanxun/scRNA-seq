#!/bin/sh
#===============================================================================
#
#          FILE:  mutect2_GATK4.sh
# 
#         USAGE:  ./mutect2_GATK4.sh 
# 
#   DESCRIPTION:  call somatic mutatioin using GATK4
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (jiangst), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  03/25/19 17:57:12 CST
#      REVISION:  ---
#===============================================================================

if [ $# -lt 1 ];then
        echo "$0 useage:"
        echo "Input the patient id:"
        echo "P101"
        exit
fi

patient=$1

#script.initialization
workdir=/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/3_Target/${patient}
tumor_bam=$workdir/${patient}TS0/bwa/${patient}TS0.sorted.markdup.BQSR.bam
normal_bam=$workdir/${patient}TSC/bwa/${patient}TSC.sorted.markdup.BQSR.bam
#reference_knownsite=/pnas/wangqf_group/database/RefDB
reference_fa=/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta
GenomeAnalysisTK4=/pnas/wangqf_group/yuxx/software/miniconda3/bin/gatk

# mutect2
rm $workdir/${patient}.mutect2*
date >>$workdir/$patient.mutect2.runlog

$GenomeAnalysisTK4 Mutect2 \
-R $reference_fa \
-I ${tumor_bam} \
-I ${normal_bam} \
-tumor ${patient}TS0 \
-normal ${patient}TSC \
-O $workdir/$patient.mutect2.vcf \
&& echo -e "mutect2 done">>$workdir/$patient.mutect2.runlog || echo -e "mutect2 crash">>$workdir/$patient.mutect2.runlog
date >> $workdir/$patient.mutect2.runlog


$GenomeAnalysisTK4 FilterMutectCalls \
-V $workdir/$patient.mutect2.vcf \
-O $workdir/$patient.mutect2.filtered.vcf \
-R $reference_fa \
&& echo -e "FilterMutectCalls done">>$workdir/$patient.mutect2.runlog || echo -e "FilterMutectCalls crash">>$workdir/$patient.mutect2.runlog
date >> $workdir/$patient.mutect2.runlog

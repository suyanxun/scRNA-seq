#!/bin/bash
#===============================================================================
#
#          FILE:  somaticsniper.sh
# 
#         USAGE:  ./somaticsniper.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  03/25/2019 09:51:48 PM CST
#      REVISION:  ---
#===============================================================================


if [ $# -lt 1 ];then
        echo "$0 useage:"
        echo "Input the patient id:"
        echo "eg: P001"
        exit
fi

patient=$1

#.script.initialization
workdir=/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/3_Target/${patient}
tumor_bam=$workdir/${patient}TS0/bwa/${patient}TS0.sorted.markdup.BQSR.bam
normal_bam=$workdir/${patient}TSC/bwa/${patient}TSC.sorted.markdup.BQSR.bam
#reference_knownsite=/pnas/wangqf_group/database/RefDB
reference_fa=/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta
#=====
Somaticsniper=/pnas/wangqf_group/database/Package/somaticsniper/somaticsniper/somatic-sniper/build/bin/bam-somaticsniper

#mutation calling using somaticsniper
date >>$workdir/${patient}.sniper.runlog
$Somaticsniper -Q 15 -G -L -F vcf \
-f $reference_fa \
$tumor_bam \
$normal_bam \
$workdir/${patient}.sniper_Q15.vcf \
&& echo -e "somaticsniper done">>$workdir/${patient}.sniper.runlog || echo -e "somaticsniper crash">>$workdir/${patient}.sniper.runlog
date >>$workdir/${patient}.sniper.runlog

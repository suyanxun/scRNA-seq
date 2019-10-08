#!/bin/bash
#===============================================================================
#
#          FILE:  annotation2furtherfilter4TargetedSeq_jiangst_.sh
# 
#         USAGE:  ./annotation2furtherfilter4TargetedSeq_jiangst_.sh 
# 
#   DESCRIPTION:  Filter the somatic mutation detected at diagnosis(mutect2,varscan,somaticsniper)
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  
#          BUGS:  ---
#         NOTES:  Need to do annotation before run the script; ARGV needed：projectname
#        AUTHOR:   (jiangst), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  04/04/2019 02:46:14 PM CST
#      REVISION:  ---
#===============================================================================


if [ $# -lt 1 ];then
	echo "$0 useage:"
	echo "Input the information in right order:"
	echo "projectname(eg:P101)"
	exit
fi

projectname=$1

#====================initiation=================
Dir=/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/3_Target/
Bam_path=/asnas/wangqf_group/jiangsht/TargetedSequencing/Result/Preprocess
reference_fa=/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta
bamreadcount=/software/biosoft/software/bam-readcount/build/bin/bam-readcount

#====================filter=====================
date +"%Y-%m-%d %H:%M:%S" > ${projectname}.annovation.filter.runlog

#select varscan_snp
perl /pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/script/target/Annovar_filter/varscan_select.pl \
-annovar_txt ${Dir}/${projectname}/${projectname}.varscan_somatic.snp.annovar.vcf.hg38_multianno.txt \
-result ${Dir}/${projectname}/${projectname}.varscan_somatic.snp.annovar.vcf.hg38_multianno.filter.txt \
&& echo -e "select varscan_snp done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "select varscan_snp crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

#select somatic sniper
perl /pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/script/target/Annovar_filter/sniper_select.pl \
-annovar_txt ${Dir}/${projectname}/${projectname}.sniper_Q15.annovar.vcf.hg38_multianno.txt \
-result ${Dir}/${projectname}/${projectname}.sniper_Q15.annovar.vcf.hg38_multianno.filter.txt \
&& echo -e "select somaticsniper done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "select somaticsniper crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

#select mutect2  
perl /pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/script/target/Annovar_filter/mutect2_select.pl \
-annovar_txt ${Dir}/${projectname}/${projectname}.mutect2.annovar.vcf.hg38_multianno.txt \
-result ${Dir}/${projectname}/${projectname}.mutect2.annovar.vcf.hg38_multianno.filter.all.txt \
-SNP ${projectname}.mutect2.annovar.vcf.hg38_multianno.filter.snp.txt \
-indel ${projectname}.mutect2.annovar.vcf.hg38_multianno.filter.indel.txt \
&& echo -e "select mutect2 done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "select mutect2 crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

#===============combine result from 3 software====================
cat ${Dir}/${projectname}/${projectname}.varscan_somatic.snp.annovar.vcf.hg38_multianno.filter.txt ${Dir}/${projectname}/${projectname}.sniper_Q15.annovar.vcf.hg38_multianno.filter.txt|sort|uniq -d > ${Dir}/${projectname}/${projectname}temp_combine_SNP.txt
cat ${Dir}/${projectname}/${projectname}temp_combine_SNP.txt ${Dir}/${projectname}/${projectname}.mutect2.annovar.vcf.hg38_multianno.filter.all.txt|sort -n -k 1.4 -k 2|uniq > ${Dir}/${projectname}/${projectname}combine_SNP.txt \
&& echo -e "combine result from 3 software done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "combine result from 3 software crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

#construct VCF
perl /pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/script/target/Annovar_filter/vcf_construct2.pl \
-txt ${Dir}/${projectname}/${projectname}combine_SNP.txt \
-vcf ${Dir}/${projectname}/${projectname}combine_SNP.vcf \
&& echo -e "construct VCF done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "construct VCF crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

#========================bamreadcount=============================
#generate hcsnp interval
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${Dir}/${projectname}/${projectname}combine_SNP.vcf > ${Dir}/${projectname}/${projectname}combine_SNP.var

##bam-readcount ; To calculate the read depth in every position
#TS0
$bamreadcount -b 15 -w 1 -l ${Dir}/${projectname}/${projectname}combine_SNP.var -f $reference_fa ${Dir}/${projectname}/${projectname}TS0/bwa/${projectname}TS0.sorted.markdup.BQSR.bam > ${Dir}/${projectname}/${projectname}TS0.combine_SNP.readcount
#TSC
$bamreadcount -b 15 -w 1 -l ${Dir}/${projectname}/${projectname}combine_SNP.var -f $reference_fa ${Dir}/${projectname}/${projectname}TSC/bwa/${projectname}TSC.sorted.markdup.BQSR.bam > ${Dir}/${projectname}/${projectname}TSC.combine_SNP.readcount

#==========calculate the VAF and combine the information==========
perl /pnas/wangqf_group/suyx/PMO/github/pipeline/scRNA-seq/bin/script/target/Annovar_filter/VAF_E0_EC.pl \
-SNP ${Dir}/${projectname}/${projectname}combine_SNP.txt \
-bamreadcount_E0 ${Dir}/${projectname}/${projectname}TS0.combine_SNP.readcount \
-bamreadcount_EC ${Dir}/${projectname}/${projectname}TSC.combine_SNP.readcount \
-output_file ${Dir}/${projectname}/${projectname}.combine_SNP.E0_EC.txt \
&& echo -e "calculate the VAF done!" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog || echo -e "calculate the VAF crash" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog

date +"%Y-%m-%d %H:%M:%S" >>${Dir}/${projectname}/${projectname}.annovation.filter.runlog


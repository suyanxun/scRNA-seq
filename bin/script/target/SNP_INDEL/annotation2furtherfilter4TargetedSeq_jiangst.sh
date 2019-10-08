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
workdir=/asnas/wangqf_group/jiangsht/TargetedSequencing/Result/SomaticMutation/Filter
Bam_path=/asnas/wangqf_group/jiangsht/TargetedSequencing/Result/Preprocess
reference_fa=/pnas/wangqf_group/database/RefDB/UCSC_hg19/Sequence/WholeGenomeFasta/genome.fa
bamreadcount=/software/biosoft/software/bam-readcount/build/bin/bam-readcount
cd $workdir

#====================filter=====================
date +"%Y-%m-%d %H:%M:%S" > ${projectname}.annovation.filter.runlog

#select varscan_snp
perl /asnas/wangqf_group/jiangsht/LDCH/Script/SomaticFilter/essential_script/varscan_select.pl \
-annovar_txt /asnas/wangqf_group/jiangsht/TargetedSequencing/Result/SomaticMutation/${projectname}.varscan_somatic.snp.vcf.annovation.hg19_multianno.txt \
-result ${projectname}.varscan_somatic.snp.vcf.annovation.hg19_multianno.filter.txt \
&& echo -e "select varscan_snp done!" >>${projectname}.annovation.filter.runlog || echo -e "select varscan_snp crash" >>${projectname}.annovation.filter.runlog

#select somatic sniper
perl /asnas/wangqf_group/jiangsht/LDCH/Script/SomaticFilter/essential_script/sniper_select.pl \
-annovar_txt /asnas/wangqf_group/jiangsht/TargetedSequencing/Result/SomaticMutation/${projectname}.sniper_Q15.vcf.annovation.hg19_multianno.txt \
-result ${projectname}.sniper_Q15.vcf.annovation.hg19_multianno.filter.txt \
&& echo -e "select somaticsniper done!" >>${projectname}.annovation.filter.runlog || echo -e "select somaticsniper crash" >>${projectname}.annovation.filter.runlog

#select mutect2  
perl /asnas/wangqf_group/jiangsht/LDCH/Script/SomaticFilter/essential_script/mutect2_select.pl \
-annovar_txt /asnas/wangqf_group/jiangsht/TargetedSequencing/Result/SomaticMutation/${projectname}.mutect2.filtered.vcf.annovation.hg19_multianno.txt \
-result ${projectname}.mutect2.filtered.vcf.annovation.hg19_multianno.filter.all.txt \
-SNP ${projectname}.mutect2.filtered.vcf.annovation.hg19_multianno.filter.SNP.txt \
-indel ${projectname}.mutect2.filtered.vcf.annovation.hg19_multianno.filter.indel.txt \
&& echo -e "select mutect2 done!" >>${projectname}.annovation.filter.runlog || echo -e "select mutect2 crash" >>${projectname}.annovation.filter.runlog

#===============combine result from 3 software====================
cat ${projectname}.varscan_somatic.snp.vcf.annovation.hg19_multianno.filter.txt ${projectname}.sniper_Q15.vcf.annovation.hg19_multianno.filter.txt|sort|uniq -d > ${projectname}temp_combine_SNP.txt
cat ${projectname}temp_combine_SNP.txt ${projectname}.mutect2.filtered.vcf.annovation.hg19_multianno.filter.all.txt|sort -n -k 1.4 -k 2|uniq > ${projectname}combine_SNP.txt \
&& echo -e "combine result from 3 software done!" >>${projectname}.annovation.filter.runlog || echo -e "combine result from 3 software crash" >>${projectname}.annovation.filter.runlog

#construct VCF
perl /asnas/wangqf_group/jiangsht/LDCH/Script/SomaticFilter/essential_script/vcf_construct2.pl \
-txt ${projectname}combine_SNP.txt \
-vcf ${projectname}combine_SNP.vcf \
&& echo -e "construct VCF done!" >>${projectname}.annovation.filter.runlog || echo -e "construct VCF crash" >>${projectname}.annovation.filter.runlog

#========================bamreadcount=============================
#generate hcsnp interval
perl -ane 'print join("\t",@F[0,1,1])."\n" unless(m/^#/)' ${projectname}combine_SNP.vcf > ${projectname}combine_SNP.var

##bam-readcount ; To calculate the read depth in every position
#TS0
$bamreadcount -b 15 -w 1 -l ${projectname}combine_SNP.var -f $reference_fa $Bam_path/${projectname}TS0.7_2.clean.sort.fixmate.dedup.indelrealgn.BQSR.bam > ${projectname}TS0.combine_SNP.readcount
#TSC
$bamreadcount -b 15 -w 1 -l ${projectname}combine_SNP.var -f $reference_fa $Bam_path/${projectname}TSC.7_2.clean.sort.fixmate.dedup.indelrealgn.BQSR.bam > ${projectname}TSC.combine_SNP.readcount

#==========calculate the VAF and combine the information==========
perl /asnas/wangqf_group/jiangsht/TargetedSequencing/Src/VAF_E0_EC.pl \
-SNP ${projectname}combine_SNP.txt \
-bamreadcount_E0 ${projectname}TS0.combine_SNP.readcount \
-bamreadcount_EC ${projectname}TSC.combine_SNP.readcount \
-output_file ${projectname}.combine_SNP.E0_EC.txt \
&& echo -e "calculate the VAF done!" >>${projectname}.annovation.filter.runlog || echo -e "calculate the VAF crash" >>${projectname}.annovation.filter.runlog

date +"%Y-%m-%d %H:%M:%S" >>${projectname}.annovation.filter.runlog


#!/bin/bash
#===============================================================================
#
#          FILE:  preprocess_1.sh
# 
#         USAGE:  ./preprocess_1.sh 
# 
#   DESCRIPTION:  Input the projectname workdir rawdatapath ppn memory:
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (jiangsht), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  06/14/2019 04:59:15 PM CST
#      REVISION:  ---
#===============================================================================


if [ $# -lt 5 ];then
        echo "$0 useage:"
        echo "Input the projectname workdir rawdatapath ppn memory:"
	echo "P101SR00"
        exit
fi


# shell script input options
projectname=$1
workdir=$2
rawdatapath=$3
ppn=$4 # the ppn
memory=$5

#output diretory
if [ ! -d $workdir/Cleanfq ]
then mkdir -p $workdir/Cleanfq
fi

if [ ! -d $workdir/Preprocess/Hg38 ]
then mkdir -p $workdir/Preprocess/Hg38
fi

#step.0.script.initialization 
rawdata_1=${projectname}_R1.fastq.gz
rawdata_2=${projectname}_R2.fastq.gz

rm -f $workdir/Preprocess/Hg38/$projectname.*
#=====
Fastqc=/asnas/wangqf_group/yuxx/software/FastQC/fastqc
Trimmomatic=/asnas/wangqf_group/yuxx/software/Trimmomatic-0.36/trimmomatic-0.36.jar
JAVA=/software/biosoft/software/jdk1.8/jdk1.8.0_45/bin/java  #java version "1.8.0_45"
BWA=/software/biosoft/software/bwa-0.7.17/bwa  #Version: 0.7.17-r1188
reference_fa=/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta
Picard=/software/biosoft/software/picard/build/libs/picard.jar
GenomeAnalysisTK=/pnas/wangqf_group/yuxx/software/gatk-4.0.11.0/gatk
GATK=/software/biosoft/software/gatk/GenomeAnalysisTK.jar 
SamTools=/software/biosoft/software/samtools1.3.1/bin/samtools  #Version: 1.9 (using htslib 1.9)
assembly=hg38
GATK_bundle=/pnas/wangqf_group/database/RefDB/GATK/resources/bundle/hg38/
TMPD=/asnas/wangqf_group/jiangsht/TargetedSequencing/Result/Temp

# trimmomatic
echo "$projectname" > $workdir/$projectname.runlog
date >> $workdir/$projectname.runlog

cd $workdir/Cleanfq
rm -f ./${projectname}.*
time $JAVA -jar $Trimmomatic PE -phred33 \
		$rawdatapath/$rawdata_1 $rawdatapath/$rawdata_2 \
	        -baseout ${projectname}.trim.fastq.gz \
		ILLUMINACLIP:${Trimmomatic%%/trimmomatic-0.36.jar}/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
		LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 -threads $ppn \
     && echo -e "** fq QC done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** fq QC crash ** @ "`date` >> $workdir/$projectname.runlog

# bowtie2.mapping
time $BWA mem -t $ppn -M -Y -R "@RG\tID:$projectname\tLB:$projectname\tSM:$projectname\tPL:ILLUMINA"  $reference_fa \
		 $workdir/Cleanfq/${projectname}.trim_1P.fastq.gz $workdir/Cleanfq/${projectname}.trim_2P.fastq.gz | $SamTools view -Sb - >$workdir/Preprocess/Hg38/${projectname}.bam \
		 && echo -e "** BWA MEM done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** BWA MEM crash ** @ "`date` >> $workdir/$projectname.runlog

# sort the bam file
#memoryP=$[$memory / $ppn]
time $SamTools sort -@ $ppn -m 2G -O bam -o $workdir/Preprocess/Hg38/${projectname}.sorted.bam $workdir/Preprocess/Hg38/${projectname}.bam \
&& echo -e "** sorted raw bamfile done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** sorted raw bamfile crash ** @ "`date` >> $workdir/$projectname.runlog

# tag the Duplicate reads
time $GenomeAnalysisTK MarkDuplicates \
		-I $workdir/Preprocess/Hg38/${projectname}.sorted.bam \
		-M $workdir/Preprocess/Hg38/${projectname}.markdup_metrics.txt \
		-O $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.bam \
     && echo "** ${projectname}.sorted.bam MarkDuplicates done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** ${projectname}.sorted.bam MarkDuplicates crash ** @ "`date` >> $workdir/$projectname.runlog
		
#build the index for ${projectname}.sorted.markdup.bam
time $SamTools index $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.bam \
&& echo -e "** ${projectname}.sorted.bam MarkDuplicates index done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** ${projectname}.sorted.bam MarkDuplicates index crash ** @ "`date` >> $workdir/$projectname.runlog

# BQSR: Base quality score recalibration
time $GenomeAnalysisTK BaseRecalibrator \
		 -R $reference_fa \
		 -I $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.bam \
		 --known-sites $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
		 --known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		 --known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
		 -O $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.recal_data.table \
     && echo -e "** ${projectname}.sorted.markdup.recal_data.table done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** ${projectname}.sorted.markdup.recal_data.table crash ** @ "`date` >> $workdir/$projectname.runlog

time $GenomeAnalysisTK ApplyBQSR \
		 --bqsr-recal-file $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.recal_data.table \
		 -R $reference_fa \
		 -I $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.bam \
		 -O $workdir/Preprocess/Hg38/${projectname}.sorted.markdup.BQSR.bam \
     && echo -e "** ${projectname}.sorted.markdup.BQSR.bam done ** @ "`date` >> $workdir/$projectname.runlog || echo -e "** ${projectname}.sorted.markdup.BQSR.bam crash ** @ "`date` >> $workdir/$projectname.runlog

date >> $workdir/$projectname.runlog

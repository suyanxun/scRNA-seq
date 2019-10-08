#!/usr/bin/bash

#this script is matched with single sample of paired fastq files that sequenced by Illumina

#the related softwares and tools
fastqc="/asnas/wangqf_group/yuxx/software/FastQC/fastqc"
trimmomatic="/asnas/wangqf_group/yuxx/software/Trimmomatic-0.36/trimmomatic-0.36.jar"
bwa="/software/biosoft/software/bwa-0.7.17/bwa"
samtools="/pnas/wangqf_group/yuxx/software/miniconda3/bin/samtools"
gatk="/pnas/wangqf_group/yuxx/software/gatk-4.0.11.0/gatk"
picard="/pnas/liuxin_group/yuxx/software/picard.jar"
GATK=""
table_annovar="/pnas/wangqf_group/database/RefDB/annovar/table_annovar.pl"
annovarfilter="/pnas/wangqf_group/database/RefDB/annovar/humandb_hg38"

#reference and files
reference="/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta"
GATK_bundle="/pnas/wangqf_group/database/RefDB/GATK/resources/bundle/hg38/"
targetRegion="/asnas/wangqf_group/yuxx/baseDataFiles/targetPanel/pediatricTargetPanel.bed"
assembly="hg38"

#it's not neccesary if there are index files for the vcf files
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_omni2.5.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/dbsnp_146.hg38.vcf
#$gatk IndexFeatureFile --feature-file $GATK_bundle/hapmap_3.3_grch38_pop_stratified_af.vcf

# shell script input options
fq1=$1
fq2=$2
RGID=$3 # Read Group, it can be replaced by Lane ID 
library=$4 # the sequence library ID
sample=$5 # the sample ID
outdir=$6 # the output path
ppn=$7 # the ppn
memory=$8

# set the diretory according the sample
outdir=${outdir}/${sample}

# get the fastq file pre-name
fq_file_name=`basename $fq1`
fq_file_name=${fq_file_name%%_1.fastq.gz}

#output diretory
if [ ! -d $outdir/cleanfq ]
then mkdir -p $outdir/cleanfq
fi

if [ ! -d $output/bwa ]
then mkdir -p $outdir/bwa
fi

if [ ! -d $output/gatk ]
then mkdir -p $outdir/gatk
fi

# use the FastQC and Trimmomatic to quality control for fastq files
cd $outdir/cleanfq
#time $fastqc -t $ppn -o ./ $fq1 $fq2 --noextract && echo "** fq fastqc done ** @ "`date` > $outdir/${sample}.checked.txt

time java -jar $trimmomatic PE -phred33 \
		 $fq1 $fq2 \
		 -baseout ${fq_file_name}.trim.fq.gz \
		 ILLUMINACLIP:${trimmomatic%%/trimmomatic-0.36.jar}/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
		 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 -threads $ppn && echo "** fq QC done ** @ "`date` >> $outdir/${sample}.checked.txt

# use the bwa to map the reference genome
time $bwa mem -t $ppn -M -Y -R "@RG\tID:$RGID\tPL:ILLUMINA\tLB:$library\tSM:$sample" $reference \
		 ${fq_file_name}.trim_1P.fq.gz ${fq_file_name}.trim_2P.fq.gz | $samtools view -Sb - > $outdir/bwa/${sample}.bam && \
		 echo "** BWA MEM done ** @ "`date` >> $outdir/${sample}.checked.txt

# sort the bam file
#memoryP=$[$memory / $ppn]
time $samtools sort -@ $ppn -m 2G -O bam -o $outdir/bwa/${sample}.sorted.bam $outdir/bwa/${sample}.bam && echo "** sorted raw bamfile done ** @ "`date` >>$outdir/${sample}.checked.txt

# tag the Duplicate reads
$gatk MarkDuplicates \
		-I $outdir/bwa/${sample}.sorted.bam \
		-M $outdir/bwa/${sample}.markdup_metrics.txt \
		-O $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done ** @ "`date` >> $outdir/${sample}.checked.txt
		
#build the index for ${sample}.sorted.markdup.bam
time $samtools index $outdir/bwa/${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates index done ** @ "`date` >> $outdir/${sample}.checked.txt

# BQSR: Base quality score recalibration
time $gatk BaseRecalibrator \
		 -R $reference \
		 -I $outdir/bwa/${sample}.sorted.markdup.bam \
		 --known-sites $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
		 --known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
		 --known-sites $GATK_bundle/dbsnp_146.hg38.vcf \
		 -O $outdir/bwa/${sample}.sorted.markdup.recal_data.table && echo "** ${sample}.sorted.markdup.recal_data.table done ** @ "`date` >> $outdir/${sample}.checked.txt

time $gatk ApplyBQSR \
		 --bqsr-recal-file $outdir/bwa/${sample}.sorted.markdup.recal_data.table \
		 -R $reference \
		 -I $outdir/bwa/${sample}.sorted.markdup.bam \
		 -O $outdir/bwa/${sample}.sorted.markdup.BQSR.bam && echo "** ${sample}.sorted.markdup.BQSR.bam done ** @ "`date` >> $outdir/${sample}.checked.txt

# Reads information assessment
# Objective: To determine the amount of available data, calculate plus measurements, calculate the proportion of duplicate data, calculate the coverage ratio of the genome, and calculate the coverage depth
echo -e "** $sample Reads information assessment starts ** @"`date` >> $outdir/${sample}.checked.txt
# Statistics for the mapping result of with removing duplicate bams
#java -jar $GATK -T FlagStat -R $reference -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam -o $outdir/bwa/${sample}.rmdup.flagstat.txt
# Statistical coverage depth information
#java -jar $GATK -T DepthOfCoverage -R $reference -I $outdir/bwa/${sample}.sorted.markdup.BQSR.bam -o $outdir/bwa/depthofcoverage
#echo -e "** $sample Reads information assessment ends @"`date` >> $outdir/${sample}.checked.txt || echo -e "** $sample Reads information assessment errors @"`date` >> $outdir/${sample}.checked.txt
#echo -e "** $sample coverage of histogram starts @"`date` >> $outdir/${sample}.checked.txt
#cat $outdir/bwa/depthofcoverage | cut -f 3 > $outdir/bwa/depthofcoverage.txt
#Rscript /pnas/wangqf_group/yuxx/scripts/Rscripts/coverageHist.R $outdir/bwa depthofcoverage.txt $sample
#echo -e "** $sample coverage of histogram ends @"`date` >> $outdir/${sample}.checked.txt || echo -e "** $sample Reads coverage of histogram errors @"`date` >> $outdir/${sample}.checked.txt


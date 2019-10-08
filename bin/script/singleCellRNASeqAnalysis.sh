#!/usr/bin/bash

#this script is used to run single sample single cell fastqs files that sequenced by Illumina

#the related softwares and tools
fastqc="/asnas/wangqf_group/yuxx/software/FastQC/fastqc"
trimmomatic="/asnas/wangqf_group/yuxx/software/Trimmomatic-0.36/trimmomatic-0.36.jar"
cutadapt="/pnas/wangqf_group/yuxx/software/miniconda3/envs/python3/bin/cutadapt"
cellranger="/pnas/wangqf_group/yuxx/software/cellranger-3.0.1/cellranger"
transcriptome="/pnas/wangqf_group/yuxx/RefenceGenome/refdata-cellranger-GRCh38-3.0.0"
# shell script input options

# Name of sample, a folder of / will be created and be the main directory
id=$1
# Sample name as specified in the sample sheet supplied to cellranger mkfastq. If have multiple libraries for the sample, it will need to run cellranger count on them individually
sample=$id
# path of Pair-end sequencing file
fastqs=/asnas/wangqf_group/suyx/Project_scRNA_seq/RawData
# Hard-trim the input R1 sequence to this length
r1_length=$2
# Hard-trim the input R2 sequence to this length
ppn=12
# memory. please set the mem(GB)
memory=60
# Expected number of recovered cells
expectCells=$3
# Output path (default ./, the result will be stored in ./)
outdir=/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/1_CellRanger_rawdata
# the remain Read2 min length after trim and cuadapt
minlen=20

#########################################################################################################################################################################################################
# set the diretory according the sample
outdir=${outdir}/${sample}

#output diretory
if [ ! -d $outdir/fastQCRes ]
then mkdir -p $outdir/fastQCRes
fi

if [ ! -d $outdir/trimRes ]
then mkdir -p $outdir/trimRes
fi

if [ ! -d $outdir/cutadapt ]
then mkdir -p $outdir/cutadapt
fi

if [ ! -d $outdir/cellrangerRes ]
then mkdir -p $outdir/cellrangerRes
fi

#######################################################################################################################################################################################################
# use the FastQC and Trimmomatic to quality control for fastq files
cd $outdir/fastQCRes
for fq2 in `ls ${fastqs}/${sample}/${sample}*_R2_001.fastq.gz`
	do
		time $fastqc -t $ppn -o ./ $fq2 --noextract && echo "** $fq2 fastqc done ** @ "`date` > $outdir/${sample}.checked.txt
	done

######################################################################################################################################################################################################
# use the Trimmomatic and cutadapt to quality control for fastq files
cd $outdir/trimRes
for fq2 in `ls ${fastqs}/${sample}/${sample}*_R2_001.fastq.gz`
	do
		name=`basename $fq2`
		name=${name%%_R2_001.fastq.gz}

		time java -jar $trimmomatic SE -phred33 \
				  $fq2 ${name}_trim_R2_001.fastq.gz \
				  ILLUMINACLIP:${trimmomatic%%/trimmomatic-0.36.jar}/adapters/TruSeq3-PE-2.fa:2:30:10 \
				  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minlen -threads $ppn && echo "** $fq2 trim done ** @ "`date` >> $outdir/${sample}.checked.txt

		time $cutadapt -a "A{100}" \
					   -o ../cutadapt/${name}_R2_001.fastq.gz ${name}_trim_R2_001.fastq.gz \
					   -m $minlen && echo "** $fq2 cutadapt done ** @ "`date` >> $outdir/${sample}.checked.txt

		time perl /pnas/wangqf_group/yuxx/scripts/Perlscripts/splitfastqFiles.pl ${fastqs}/${sample}/${name}_R1_001.fastq.gz $outdir/cutadapt/${name}_R2_001.fastq.gz $outdir/cutadapt/${name}_R1_001.fastq.gz
		echo "** ${name}_R1_001.fastq.gz correct done ** @ "`date` >> $outdir/${sample}.checked.txt

	done

#####################################################################################################################################################################################################

fastqs=$outdir/cutadapt
cd $outdir/cellrangerRes
$cellranger count --id=$id \
				  --sample=$sample \
		  		  --transcriptome=$transcriptome \
			      --fastqs=$fastqs \
				  --r1-length=$r1_length \
			      --localcores=$ppn \
				  --localmem=$memory 
				  --expect-cells=$expectCells \

echo "** $sample cellranger count done ** @ "`date` >> $outdir/${sample}.checked.txt

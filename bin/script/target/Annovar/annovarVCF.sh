#!/usr/bin/bash

#this script is used for annovar varients

table_annovar="/pnas/wangqf_group/database/RefDB/annovar/table_annovar.pl"
annovarfilter="/pnas/wangqf_group/database/RefDB/annovar/humandb_hg38"

#reference and files
reference="/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta"
assembly="hg38"

sample=$1
vcfFile=$2
outdir=$3 # the output path
ppn=$4 # the ppn

# annotation for SNP or indel use annovar
perl $table_annovar $vcfFile $annovarfilter \
		 -buildver $assembly \
		 -out $outdir/${sample}.annovar.vcf \
		 -thread $ppn \
		 -remove -protocol refGene,cytoBand,avsnp150,dbnsfp35c,cosmic87_coding,cosmic87_noncoding,clinvar_20180603,ALL.sites.2015_08,AFR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,EUR.sites.2015_08,SAS.sites.2015_08,dbscsnv11,esp6500siv2_all,exac03,exac03nonpsych,genomicSuperDups,gnomad_genome,intervar_20180118,kaviar_20150923,nci60 \
		 -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
		 -vcfinput
echo -e "** the $sample has finished the annotation ** @ "`date` >> $outdir/$sample.checked.txt || cat -e "** the $sample annotated errors ** @ "`date` >> $outdir/$sample.checked.txt


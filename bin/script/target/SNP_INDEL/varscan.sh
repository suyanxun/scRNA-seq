#!/bin/bash
#===============================================================================
#
#          FILE:  varscan.sh
# 
#         USAGE:  ./varscan.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (jiangst), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  03/25/2019 08:17:52 PM CST
#      REVISION:  ---
#===============================================================================

if [ $# -lt 1 ];then
        echo "$0 useage:"
        echo "Input the patient id:"
        echo "P101"
        exit
fi

patient=$1

#step.0.script.initialization
workdir=/asnas/wangqf_group/suyx/Project_scRNA_seq/Analysis/Result/3_Target/${patient}
tumor_bam=$workdir/${patient}TS0/bwa/${patient}TS0.sorted.markdup.BQSR.bam
normal_bam=$workdir/${patient}TSC/bwa/${patient}TSC.sorted.markdup.BQSR.bam
#reference_knownsite=/pnas/wangqf_group/database/RefDB
reference_fa=/pnas/wangqf_group/database/RefDB/GATK/reference/hg38/Homo_sapiens_assembly38.fasta
interval_file=/asnas/wangqf_group/jiangsht/TargetedSequencing/Src/targetedRegionFromRefgene.sort.withoutGene.bed

#=====
JAVA=/software/biosoft/software/jdk1.8/jdk1.8.0_45/bin/java
VarScan=/software/biosoft/software/varscan/VarScan.v2.4.3.jar
SamTools=/pnas/wangqf_group/database/Package/Samtools/samtools-1.8/samtools  #Version: 1.9 (using htslib 1.9)

#step.11 varscan_somatic
#tumor mpileup
date > $workdir/$patient.varscan.runlog
$SamTools  mpileup -A -f $reference_fa -l $interval_file  ${tumor_bam} > $workdir/${patient}TS0.sorted.markdup.BQSR.exon.mpileup   \
&& perl -lane 'BEGIN{$k1=0;$k10=0;$k20=0;$k30=0;$k40=0;$k50=0;$k60=0;$k70=0;$k80=0;$k90=0;$k100=0;$k110=0;$k120=0;$k130=0;$k140=0;$k150=0;$k160=0;$k170=0;$k180=0;$k190=0;$k200=0;$k210=0;$k220=0;$k230=0;$k240=0;$k250=0;$k260=0;$k270=0;$k280=0;$k290=0;$k300=0;$k310=0;$k320=0;$k330=0;$k340=0;$k350=0;$k360=0;$k370=0;$k380=0;$k390=0;$k400=0;$k410=0;$k420=0;$k430=0;$k440=0;$k450=0;$k460=0;$k470=0;$k480=0;$k490=0;$k500=0;}$a=$F[3];if($a >=1){$k1 ++};if($a >=10){$k10++};if($a >=20){$k20++};if($a >=30){$k30++};if($a >=40){$k40++};if($a >=50){$k50++};if($a >=60){$k60++};if($a >=70){$k70++};if($a >=80){$k80++};if($a >=90){$k90++};if($a >=100){$k100++};if($a >=110){$k110++};if($a >=120){$k120++};if($a >=130){$k130++};if($a >=140){$k140++};if($a >=150){$k150++};if($a >=160){$k160++};if($a >=170){$k170++};if($a >=180){$k180++};if($a >=190){$k190++};if($a >=200){$k200++};if($a >=210){$k210++};if($a >=220){$k220++};if($a >=230){$k230++};if($a >=240){$k240++};if($a >=250){$k250++};if($a >=260){$k260++};if($a >=270){$k270++};if($a >=280){$k280++};if($a >=290){$k290++};if($a >=300){$k300++};if($a >=310){$k310++};if($a >=320){$k320++};if($a >=330){$k330++};if($a >=340){$k340++};if($a >=350){$k350++};if($a >=360){$k360++};if($a >=370){$k370++};if($a >=380){$k380++};if($a >=390){$k390++};if($a >=400){$k400++};if($a >=410){$k410++};if($a >=420){$k420++};if($a >=430){$k430++};if($a >=440){$k440++};if($a >=450){$k450++};if($a >=460){$k460++};if($a >=470){$k470++};if($a >=480){$k480++};if($a >=490){$k490++};if($a >=500){$k500++};END{print(join("\t",">=1X",">=10X",">=20X",">=30X",">=40X",">=50X",">=60X",">=70X",">=80X",">=90X",">=100X",">=200X",">=210X",">=220X",">=230X",">=240X",">=250X",">=260X",">=270X",">=280X",">=290X",">=300X",">=310X",">=320X",">=330X",">=340X",">=350X",">=360X",">=370X",">=380X",">=390X",">=400X",">=410X",">=420X",">=430X",">=440X",">=450X",">=460X",">=470X",">=480X",">=490X",">=500X"));print(join("\t",$k1,$k10,$k20,$k30,$k40,$k50,$k60,$k70,$k80,$k90,$k100,$k110,$k120,$k130,$k140,$k150,$k160,$k170,$k180,$k190,$k200,$k210,$k220,$k230,$k240,$k250,$k260,$k270,$k280,$k290,$k300,$k310,$k320,$k330,$k340,$k350,$k360,$k370,$k380,$k390,$k400,$k410,$k420,$k430,$k440,$k450,$k460,$k470,$k480,$k490,$k500))}'  \
$workdir/${patient}TS0.sorted.markdup.BQSR.exon.mpileup > $workdir/${patient}TS0.sorted.markdup.BQSR.exon.DPdistribution \
&& echo -e "tumor mpileup targeted region DPdistribution done!" >>$workdir/$patient.varscan.runlog || echo -e "tumor mpileup targeted region DPdistribution crash" >>$workdir/$patient.varscan.runlog 

#normal mpileup
date >> $workdir/$patient.varscan.runlog
$SamTools  mpileup -A -f $reference_fa -l $interval_file  ${normal_bam} > $workdir/${patient}TSC.sorted.markdup.BQSR.exon.mpileup   \
&& perl -lane 'BEGIN{$k1=0;$k10=0;$k20=0;$k30=0;$k40=0;$k50=0;$k60=0;$k70=0;$k80=0;$k90=0;$k100=0;$k110=0;$k120=0;$k130=0;$k140=0;$k150=0;$k160=0;$k170=0;$k180=0;$k190=0;$k200=0;$k210=0;$k220=0;$k230=0;$k240=0;$k250=0;$k260=0;$k270=0;$k280=0;$k290=0;$k300=0;$k310=0;$k320=0;$k330=0;$k340=0;$k350=0;$k360=0;$k370=0;$k380=0;$k390=0;$k400=0;$k410=0;$k420=0;$k430=0;$k440=0;$k450=0;$k460=0;$k470=0;$k480=0;$k490=0;$k500=0;}$a=$F[3];if($a >=1){$k1 ++};if($a >=10){$k10++};if($a >=20){$k20++};if($a >=30){$k30++};if($a >=40){$k40++};if($a >=50){$k50++};if($a >=60){$k60++};if($a >=70){$k70++};if($a >=80){$k80++};if($a >=90){$k90++};if($a >=100){$k100++};if($a >=110){$k110++};if($a >=120){$k120++};if($a >=130){$k130++};if($a >=140){$k140++};if($a >=150){$k150++};if($a >=160){$k160++};if($a >=170){$k170++};if($a >=180){$k180++};if($a >=190){$k190++};if($a >=200){$k200++};if($a >=210){$k210++};if($a >=220){$k220++};if($a >=230){$k230++};if($a >=240){$k240++};if($a >=250){$k250++};if($a >=260){$k260++};if($a >=270){$k270++};if($a >=280){$k280++};if($a >=290){$k290++};if($a >=300){$k300++};if($a >=310){$k310++};if($a >=320){$k320++};if($a >=330){$k330++};if($a >=340){$k340++};if($a >=350){$k350++};if($a >=360){$k360++};if($a >=370){$k370++};if($a >=380){$k380++};if($a >=390){$k390++};if($a >=400){$k400++};if($a >=410){$k410++};if($a >=420){$k420++};if($a >=430){$k430++};if($a >=440){$k440++};if($a >=450){$k450++};if($a >=460){$k460++};if($a >=470){$k470++};if($a >=480){$k480++};if($a >=490){$k490++};if($a >=500){$k500++};END{print(join("\t",">=1X",">=10X",">=20X",">=30X",">=40X",">=50X",">=60X",">=70X",">=80X",">=90X",">=100X",">=200X",">=210X",">=220X",">=230X",">=240X",">=250X",">=260X",">=270X",">=280X",">=290X",">=300X",">=310X",">=320X",">=330X",">=340X",">=350X",">=360X",">=370X",">=380X",">=390X",">=400X",">=410X",">=420X",">=430X",">=440X",">=450X",">=460X",">=470X",">=480X",">=490X",">=500X"));print(join("\t",$k1,$k10,$k20,$k30,$k40,$k50,$k60,$k70,$k80,$k90,$k100,$k110,$k120,$k130,$k140,$k150,$k160,$k170,$k180,$k190,$k200,$k210,$k220,$k230,$k240,$k250,$k260,$k270,$k280,$k290,$k300,$k310,$k320,$k330,$k340,$k350,$k360,$k370,$k380,$k390,$k400,$k410,$k420,$k430,$k440,$k450,$k460,$k470,$k480,$k490,$k500))}'  \
$workdir/${patient}TSC.sorted.markdup.BQSR.exon.mpileup > $workdir/${patient}TSC.sorted.markdup.BQSR.exon.DPdistribution \
&& echo -e "normal mpileup targeted region DPdistribution done!" >>$workdir/$patient.varscan.runlog || echo -e "normal mpileup targeted region DPdistribution crash" >>$workdir/$patient.varscan.runlog 

#varscan
date >> $workdir/$patient.varscan.runlog
$JAVA -jar $VarScan somatic \
$workdir/${patient}TSC.sorted.markdup.BQSR.exon.mpileup \
$workdir/${patient}TS0.sorted.markdup.BQSR.exon.mpileup \
--output-snp $workdir/${patient}.varscan_somatic.snp \
--output-indel $workdir/${patient}.varscan_somatic.indel \
--strand-filter 1 \
--output-vcf 1 \
--min-var-freq 0.01 \
&& echo -e "varscan_somatic done">>$workdir/$patient.varscan.runlog|| echo -e "varscan_somatic crash" >> $workdir/$patient.varscan.runlog
date >> $workdir/$patient.varscan.runlog

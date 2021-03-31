#!/bin/bash
#SBATCH --job-name=re-run_DD0001_WES         # Job name (DD0001_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=8           # Number of CPU cores per task (4)
#SBATCH --mem=70G                   # Job memory limit (10 GB)
#SBATCH --time=100:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=DD0001_WES.%j.out    # Standard output log
#SBATCH --error=DD0001_WES.%j.err     # Standard error log

thread=8
Bioproject='PRJEB12081'
Normal_Run='ERR1681499-ERR1681454'
Tumor_Run='ERR1681500-ERR1681455'
SampleName='DD0001'
results='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/results/'${SampleName} #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/data/'${SampleName}'/'${Normal_Run}
dataT='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/data/'${SampleName}'/'${Tumor_Run}
MuTect_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/'${SampleName}''
MuTect2_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect2/'${SampleName}''
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
# Germline_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/Germline/'${SampleName}''
# DepthOfCoverage='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/DepthOfCoverage/'${SampleName}''
strelka_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/strelka/'${SampleName}''
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
scripts='/home/kh31516/kh31516_Lab_Share_script'


module load SAMtools/1.9-GCC-8.3.0
module load SRA-Toolkit/2.9.6-1-centos_linux64
module load picard/2.21.6-Java-11
ml GATK/3.8-1-Java-1.8.0_144
module load Perl/5.26.1-GCCcore-6.4.0
# ml BWA/0.7.17-GCC-8.3.0

mkdir -p $MuTect_out
mkdir -p $MuTect2_out
mkdir -p $strelka_out
# ## merge_bam ##
# samtools merge ${results}/${Normal_Run}.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-1/ERR1681499.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-2/ERR1681454.bam
# samtools merge ${results}/${Tumor_Run}.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-1/ERR1681500.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-2/ERR1681455.bam  


####### Picard #######
# get header information
samtools view -H ${results}/${Normal_Run}.bam > $results/header_N
samtools view -H ${results}/${Tumor_Run}.bam > $results/header_T

# exclude unmapped reads based on FLAG
samtools view $results/${Normal_Run}.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $results/${Normal_Run}-cleaned.sam
cat $results/header_N $results/${Normal_Run}-cleaned.sam > $results/fooN
mv $results/fooN $results/${Normal_Run}-cleaned.sam
samtools view $results/${Tumor_Run}.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $results/${Tumor_Run}-cleaned.sam
cat $results/header_T $results/${Tumor_Run}-cleaned.sam > $results/fooT
mv $results/fooT $results/${Tumor_Run}-cleaned.sam

####### Picard #######
# picard sort
cd ${results}
module load picard/2.21.6-Java-11

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Normal_Run}-cleaned.sam \
O=${results}/${Normal_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Normal_Run}_normal

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Tumor_Run}-cleaned.sam \
O=${results}/${Tumor_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Tumor_Run}_tumor

# picard MarkDuplicates
java -jar $EBROOTPICARD/picard.jar  MarkDuplicates I=${results}/${Normal_Run}_rg_added_sorted.bam \
O=${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Normal_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar  MarkDuplicates I=${results}/${Tumor_Run}_rg_added_sorted.bam \
O=${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Tumor_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

# ####### Remove unneeded files #######
# rm $results/${Normal_Run}.sam $results/ERR1681500.sam 
# rm $results/header_N $results/header_T
# rm $results/${Normal_Run}_rg_added_sorted.bam $results/${Tumor_Run}_rg_added_sorted.bam
# rm $results/${Normal_Run}-cleaned.sam $results/${Tumor_Run}-cleaned.sam

ml GATK/3.8-1-Java-1.8.0_144
####### GATK Realign #######
# Generating interval file for sort.bam

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

# ######## Mutect #######
# # Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
# cd $MuTect_out
# module load MuTect/1.1.7-Java-1.7.0_80

# time java -jar $EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect \
# --reference_sequence ${reference}/canFam3.fa \
# --dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
# --defaultBaseQualities 30 \
# --intervals ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
# --input_file:normal ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
# --input_file:tumor ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
# --out ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
# --coverage_file ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
# --vcf ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf

# ####### Mutect2 #######
# #time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
# #-R $reference/canFam3.fa \
# #-I:tumor $results/ERR1681500_rg_added_sorted_dedupped_removed.realigned.bam \
# #-I:normal $results/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned \
# #--dbsnp $reference/dbsnp_9615.vcf \
# #-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
# #-o $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf



# ####### Annovar #######
# # Extract PASS records from vcf
# module load Perl/5.26.1-GCCcore-6.4.0

# awk '$7 == "PASS" {print $0}' ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# # annovar input preparation
# perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# # annovar annotate
# perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# # add gene name
# python2 $scripts/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



# ### Sanger 5 steps filtering ###
# # 5 Steps filtering

# ## With 5 steps filtering 

# ml Anaconda3/2020.02
# grep -w KEEP ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > ${MuTect_out}/${SampleName}_PASS.stat

# python $scripts/Filter_MutectStat_5steps.py \
# ${MuTect_out}/${SampleName}_PASS.stat \
# ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
# ${MuTect_out}/${SampleName}_vaf_before.txt \
# ${MuTect_out}/${SampleName}_vaf_after.txt \
# ${MuTect_out}/${SampleName}_whyout.txt \
# $SampleName 

# ## 5 steps filtering created XXX_filterMut file 

# # annovar input preparation
# perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# # annovar annotate
# perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# # add gene name
# python2 $scripts/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



####### Mutect2 #######

cd ${MuTect2_out}

module load SAMtools/1.9-GCC-8.3.0
module load GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

Normal_sample=$(samtools view -H ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--callable-depth 8 \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${MuTect2_out}/${SampleName}-f1r2.tar.gz \
--dont-use-soft-clipped-bases true \
-L ${MuTect2_source}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${MuTect2_out}/${SampleName}_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/${SampleName}-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
o=-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/${SampleName}_MuTect2_GATK4_noDBSNP.vcf \
O=${MuTect2_out}/filtered-${SampleName}_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -Xmx32g -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-${SampleName}_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf


## with PON
java -Xmx32g -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf


### Mutect2 5 steps filtering

cd ${MuTect2_out}

ml Anaconda3/2020.02
#source activate py35

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/${SampleName}_VAF_Before.txt \
${MuTect2_out}/${SampleName}_VAF_After.txt \
${MuTect2_out}/${SampleName}_whyout.txt \
${SampleName}


## Strelka
cd ${results}
module load SAMtools/1.9-GCC-8.3.0
#samtools index ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam
#samtools index ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

${script}/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--tumorBam ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--referenceFasta ${reference}/canFam3.fa \
--runDir $strelka_out

# your result will be generated in demo_somatic/results/variants
$strelka_out/runWorkflow.py -m local -j 20


## limit strelka result into CDS region

python2 $script/Limit_vcf_to_CDS.py $strelka_out/results/variants/somatic.indels.vcf.gz $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list


### Annovar for strelka
cd $strelka_out/results/

####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS > $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS

# annovar input preparation
perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS > \
$strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput

# annovar annotate
perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput \
$reference/annovar_CanFam3.1.99.gtf

# add gene name
python2 $script/Add_GeneName_N_Signature.py $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

# ####### remove unneeded files #######
# rm $results/*_sorted_dedupped_removed.bam
# rm $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
# rm $results/*.bai
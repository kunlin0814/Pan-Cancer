#!/bin/bash
#PBS -N OSA1
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=60gb
#PBS -q batch

results='/scratch/kh31516/Pan_cancer/Osteo/results/PRJEB7540/OSA1' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Pan_cancer/Osteo/data/PRJEB7540/SAMEA2706864'
dataT='/scratch/kh31516/Pan_cancer/Osteo/data/PRJEB7540/SAMEA2706865'
dataT_1='/scratch/kh31516/Pan_cancer/Osteo/data/PRJEB7540/SAMEA2706866'
tumor_combine='/scratch/kh31516/Pan_cancer/Osteo/data/PRJEB7540/SAMEA2706865-SAMEA2706866'
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/Lab_shared_PanCancer/script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Mutect/OSA1'
MuTect2_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Mutect2/OSA1'
Germline_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Germline/OSA1'
DepthOfCoverage='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/DepthOfCoverage/OSA1'
sraRunTable='/scratch/kh31516/Pan_cancer/Osteo/source/PRJEB7540_OM.txt'
download_sra_sample_scripts='/work/szlab/kh31516_Lab_Share_script'



module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0

mkdir -p $results
mkdir -p $dataN
mkdir -p $dataT
mkdir -p $dataT_1
mkdir -p ${tumor_combine}
mkdir -p ${MuTect_out}
mkdir -p ${MuTect2_out}
mkdir -p ${Germline_out}
mkdir -p ${DepthOfCoverage}

####### Download #######


cd $dataN
#fastq-dump --split-files --gzip SRR7780735
$download_sra_sample_scripts/download_sra_sample.sh SAMEA2706864 $sraRunTable
cd $dataT
#fastq-dump --split-files --gzip SRR7780736
$download_sra_sample_scripts/download_sra_sample.sh SAMEA2706865 $sraRunTable

cd $dataT_1
#fastq-dump --split-files --gzip SRR7780736
$download_sra_sample_scripts/download_sra_sample.sh SAMEA2706866 $sraRunTable

zcat $dataT/SAMEA2706865_1.fastq.gz $dataT_1/SAMEA2706866_1.fastq.gz > ${tumor_combine}/SAMEA2706865-SAMEA2706866_1.fastq.gz
zcat $dataT/SAMEA2706865_2.fastq.gz $dataT_1/SAMEA2706866_2.fastq.gz > ${tumor_combine}/SAMEA2706865-SAMEA2706866_2.fastq.gz



cd ${results}
module load BWA/0.7.17-foss-2016b

# Tumor
time bwa aln ${reference}/canFam3.fa ${tumor_combine}/SAMEA2706865-SAMEA2706866_1.fastq.gz > ${results}/SAMEA2706865-SAMEA2706866_1.sai
time bwa aln ${reference}/canFam3.fa ${tumor_combine}/SAMEA2706865-SAMEA2706866_2.fastq.gz > ${results}/SAMEA2706865-SAMEA2706866_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/SAMEA2706865-SAMEA2706866_1.sai \
${results}/SAMEA2706865-SAMEA2706866_2.sai \
${tumor_combine}/SAMEA2706865-SAMEA2706866_1.fastq.gz \
${tumor_combine}/SAMEA2706865-SAMEA2706866_2.fastq.gz \
> ${results}/SAMEA2706865-SAMEA2706866.sam

# Normal
time bwa aln ${reference}/canFam3.fa $dataN/SAMEA2706864_1.fastq.gz > ${results}/SAMEA2706864_1.sai
time bwa aln ${reference}/canFam3.fa $dataN/SAMEA2706864_2.fastq.gz > ${results}/SAMEA2706864_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/SAMEA2706864_1.sai \
${results}/SAMEA2706864_2.sai \
$dataN/SAMEA2706864_1.fastq.gz \
$dataN/SAMEA2706864_2.fastq.gz \
> ${results}/SAMEA2706864.sam

####### Convert sam to bam file #######

module load SAMtools/1.9-foss-2016b

samtools view -bS ${results}/SAMEA2706865-SAMEA2706866.sam > ${results}/SAMEA2706865-SAMEA2706866.bam
samtools view -bS ${results}/SAMEA2706864.sam > ${results}/SAMEA2706864.bam

# get header information
#grep '@SQ\|@PG' ${results}/SAMEA2706864.sam > ${results}/header_N
#grep '@SQ\|@PG' ${results}/SAMEA2706865-SAMEA2706866.sam > ${results}/header_T
samtools view -H ${results}/SAMEA2706864.bam > ${results}/header_N
samtools view -H ${results}/SAMEA2706865-SAMEA2706866.bam > ${results}/header_T

# exclude unmapped reads based on FLAG
cat ${results}/SAMEA2706864.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/SAMEA2706864-cleaned.sam
cat ${results}/header_N ${results}/SAMEA2706864-cleaned.sam > ${results}/fooN
mv ${results}/fooN ${results}/SAMEA2706864-cleaned.sam

cat ${results}/SAMEA2706865-SAMEA2706866.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/SAMEA2706865-SAMEA2706866-cleaned.sam
cat ${results}/header_T ${results}/SAMEA2706865-SAMEA2706866-cleaned.sam > ${results}/fooT
mv ${results}/fooT ${results}/SAMEA2706865-SAMEA2706866-cleaned.sam

####### Picard #######
# picard sort
cd ${results}
module load picard/2.16.0-Java-1.8.0_144

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${results}/SAMEA2706864-cleaned.sam \
O=${results}/SAMEA2706864_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=SAMEA2706864_normal

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${results}/SAMEA2706865-SAMEA2706866-cleaned.sam \
O=${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=SAMEA2706865-SAMEA2706866_tumor

# picard MarkDuplicates
java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${results}/SAMEA2706864_rg_added_sorted.bam \
O=${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/SAMEA2706864-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted.bam \
O=${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/SAMEA2706865-SAMEA2706866-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### GATK Realign #######
# Generating interval file for sort.bam
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 

cd ${MuTect_out}
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${MuTect_out}/OSA1_PASS.stat

python $script/Filter_MutectStat_5steps.py ${MuTect_out}/OSA1_PASS.stat ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



############################## Germline mutation preparation Start ######################################
# GATK 3 
module load GATK/3.8-1-Java-1.8.0_144

cd $Germline_out

# Variant calling
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

cd ${Germline_out}
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation
module load Perl/5.26.1-GCCcore-6.4.0

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output

cd ${Germline_out}
python $script/Add_GeneName_N_Signature.py ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python $script/Add_GeneName_N_Signature.py ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


#### GATK callable ####
cd ${Germline_out}
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/SAMEA2706864_table.txt \
-o ${Germline_out}/SAMEA2706864_callable_status.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/SAMEA2706865-SAMEA2706866_table.txt \
-o ${Germline_out}/SAMEA2706865-SAMEA2706866_callable_status.bed 

############################## Germline mutation preparation End ######################################

#### GATK DepthofCoverage ####
cd ${DepthOfCoverage}
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/SAMEA2706864_DepthofCoverage_CDS.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/SAMEA2706865-SAMEA2706866_DepthofCoverage_CDS.bed 



####### Mutect2 #######

cd ${MuTect2_out}
module load SAMtools/1.9-foss-2016b
module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

Normal_sample=$(samtools view -H ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${results}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--callable-depth 8 \
--f1r2-tar-gz ${MuTect2_out}/OSA1-f1r2.tar.gz \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
-L ${MuTect2_source}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${MuTect2_out}/OSA1_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/OSA1-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/OSA1_MuTect2_GATK4_noDBSNP.vcf \
-O ${MuTect2_out}/filtered-OSA1_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -cp  $script/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-OSA1_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-OSA1_MuTect2_GATK4.vcf


## with PON
java -cp  $script/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-OSA1_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-OSA1_MuTect2_GATK4.vcf


### Mutect2 5 steps filtering

cd ${MuTect2_out}

module load Anaconda3/2018.12
source activate py35

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-OSA1_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-OSA1_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-OSA1_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-OSA1_MuTect2_GATK4.vcf \
${MuTect2_out}/OSA1_VAF_Before.txt \
${MuTect2_out}/OSA1_VAF_After.txt



: "
####### remove unneeded files #######
rm ${results}/*_sorted_dedupped_removed.bam
rm ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS ${MuTect_out}/OSA1_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm ${results}/*.bai
rm ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/SAMEA2706864_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/SAMEA2706865-SAMEA2706866_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${results}/SAMEA2706864.sam ${results}/SAMEA2706865-SAMEA2706866.sam 
rm ${results}/header_N ${results}/header_T
rm ${results}/SAMEA2706864_rg_added_sorted.bam ${results}/SAMEA2706865-SAMEA2706866_rg_added_sorted.bam
rm ${results}/SAMEA2706864-cleaned.sam ${results}/SAMEA2706865-SAMEA2706866-cleaned.sam
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm ${results}/*.sai
rm ${MuTect2_out}/DbSNP_filtered-OSA1_MuTect2_GATK4.vcf
"
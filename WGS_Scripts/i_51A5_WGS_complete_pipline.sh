#!/bin/bash
#PBS -N i_51A5
#PBS -l walltime=650:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=90gb
#PBS -q batch


################################################################################################################################
###### complete pipeline (picard, GATK, annova, Mutect, Mutect2,Strelka, Germline mutation, Depth of coverage, Callablebases) ##########
##### WGS analysis needs to change the interval region #####

Normal_Run='SRR10351647'
Tumor_Run='SRR10351648'
results='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/results/i_51A5' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/i_51A5/'${Normal_Run}
dataT='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/i_51A5/'${Tumor_Run}
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect/i_51A5'
MuTect2_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect2/i_51A5'
Germline_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Germline/i_51A5'
DepthOfCoverage='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/DepthOfCoverage/i_51A5'
strelka_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/strelka/'
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


####### Download #######
module load SRA-Toolkit/2.9.1-centos_linux64

mkdir -p $dataT
mkdir -p $dataN
mkdir -p ${results}
mkdir -p ${MuTect_out}
mkdir -p ${MuTect2_out}
mkdir -p ${Germline_out}
mkdir -p ${DepthOfCoverage}
mkdir -p $strelka_out

cd $dataN

fastq-dump --split-files --gzip ${Normal_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Normal_Run} $sraRunTable

cd $dataT

fastq-dump --split-files --gzip ${Tumor_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Tumor_Run} $sraRunTable

####### BWA Mapping #######
cd ${results}
module load BWA/0.7.17-foss-2016b

# Tumor
time bwa aln -t 4 ${reference}/canFam3.fa $dataT/${Tumor_Run}_1.fastq.gz > ${results}/${Tumor_Run}_1.sai
time bwa aln -t 4 ${reference}/canFam3.fa $dataT/${Tumor_Run}_2.fastq.gz > ${results}/${Tumor_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Tumor_Run}_1.sai \
${results}/${Tumor_Run}_2.sai \
$dataT/${Tumor_Run}_1.fastq.gz \
$dataT/${Tumor_Run}_2.fastq.gz \
> ${results}/${Tumor_Run}.sam

# Normal
time bwa aln -t 4 ${reference}/canFam3.fa $dataN/${Normal_Run}_1.fastq.gz > ${results}/${Normal_Run}_1.sai
time bwa aln -t 4 ${reference}/canFam3.fa $dataN/${Normal_Run}_2.fastq.gz > ${results}/${Normal_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Normal_Run}_1.sai \
${results}/${Normal_Run}_2.sai \
$dataN/${Normal_Run}_1.fastq.gz \
$dataN/${Normal_Run}_2.fastq.gz \
> ${results}/${Normal_Run}.sam


####### Convert sam to bam file #######

cd ${results}
module load SAMtools/1.9-foss-2016b

samtools view -bS ${results}/${Tumor_Run}.sam > ${results}/${Tumor_Run}.bam
samtools view -bS ${results}/${Normal_Run}.sam > ${results}/${Normal_Run}.bam

# get header information
#grep '@SQ\|@PG' ${results}/${Normal_Run}.sam > ${results}/header_N
#grep '@SQ\|@PG' ${results}/${Tumor_Run}.sam > ${results}/header_T
samtools view -H ${results}/${Normal_Run}.bam > ${results}/header_N
samtools view -H ${results}/${Tumor_Run}.bam > ${results}/header_T

# exclude unmapped reads based on FLAG
cat ${results}/${Normal_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Normal_Run}-cleaned.sam
cat ${results}/header_N ${results}/${Normal_Run}-cleaned.sam > ${results}/fooN
mv ${results}/fooN ${results}/${Normal_Run}-cleaned.sam

cat ${results}/${Tumor_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Tumor_Run}-cleaned.sam
cat ${results}/header_T ${results}/${Tumor_Run}-cleaned.sam > ${results}/fooT
mv ${results}/fooT ${results}/${Tumor_Run}-cleaned.sam

####### Picard #######
# picard sort
cd ${results}
module load picard/2.16.0-Java-1.8.0_144

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${results}/${Normal_Run}-cleaned.sam \
O=${results}/${Normal_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Normal_Run}_normal

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${results}/${Tumor_Run}-cleaned.sam \
O=${results}/${Tumor_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Tumor_Run}_tumor

# picard MarkDuplicates
java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${results}/${Normal_Run}_rg_added_sorted.bam \
O=${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Normal_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${results}/${Tumor_Run}_rg_added_sorted.bam \
O=${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Tumor_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### GATK Realign #######
# Generating interval file for sort.bam
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
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

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 

cd ${MuTect_out}
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.interval_list \
--input_file:normal ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${MuTect_out}/i_51A5_PASS.stat

python $script/Filter_MutectStat_5steps.py ${MuTect_out}/i_51A5_PASS.stat ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


############################## Germline mutation preparation Start ######################################
# GATK 3 
module load GATK/3.8-1-Java-1.8.0_144

cd $Germline_out

# Variant calling
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

cd ${Germline_out}
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation
module load Perl/5.26.1-GCCcore-6.4.0

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output

cd ${Germline_out}
python $script/Add_GeneName_N_Signature.py ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python $script/Add_GeneName_N_Signature.py ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


#### GATK callable ####
cd ${Germline_out}
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Normal_Run}_table.txt \
-o ${Germline_out}/${Normal_Run}_callable_status.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Tumor_Run}_table.txt \
-o ${Germline_out}/${Tumor_Run}_callable_status.bed 

############################## Germline mutation preparation End ######################################


#### GATK DepthofCoverage ####
cd ${DepthOfCoverage}
module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Normal_Run}_DepthofCoverage_CDS.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Tumor_Run}_DepthofCoverage_CDS.bed 


####### Mutect2 #######

cd ${MuTect2_out}
module load SAMtools/1.9-foss-2016b
module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

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
--f1r2-tar-gz ${MuTect2_out}/i_51A5-f1r2.tar.gz \
-L ${MuTect2_source}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${MuTect2_out}/i_51A5_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/i_51A5-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/i_51A5_MuTect2_GATK4_noDBSNP.vcf \
-O ${MuTect2_out}/filtered-i_51A5_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-i_51A5_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf


## with PON
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf


### Mutect2 5 steps filtering

cd ${MuTect2_out}

module load Anaconda3/2018.12
source activate py35

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf \
${MuTect2_out}/i_51A5_VAF_Before.txt \
${MuTect2_out}/i_51A5_VAF_After.txt


## Strelka
cd ${results}
module load SAMtools/1.9-foss-2016b 
#samtools index ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam
#samtools index ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

${script}/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--tumorBam ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--referenceFasta ${reference}/canFam3.fa \
--runDir $strelka_out/i_51A5

# your result will be generated in demo_somatic/results/variants
$strelka_out/i_51A5/runWorkflow.py -m local -j 20


## limit strelka result into CDS region

python $script/Limit_vcf_to_CDS.py $strelka_out/results/variants/somatic.indels.vcf.gz ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list




: "
####### remove unneeded files #######
rm ${results}/*_sorted_dedupped_removed.bam
rm ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS ${MuTect_out}/i_51A5_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm ${results}/*.bai
rm ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${results}/${Normal_Run}.sam ${results}/${Tumor_Run}.sam 
rm ${results}/header_N ${results}/header_T
rm ${results}/${Normal_Run}_rg_added_sorted.bam ${results}/${Tumor_Run}_rg_added_sorted.bam
rm ${results}/${Normal_Run}-cleaned.sam ${results}/${Tumor_Run}-cleaned.sam
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm ${results}/*.sai
rm ${MuTect2_out}/DbSNP_filtered-i_51A5_MuTect2_GATK4.vcf
rm ${MuTect2_out}/read-orientation-model.tar.gz
rm ${MuTect2_out}/i_51A5-f1r2.tar.gz
"

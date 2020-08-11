#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=550:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=60gb
#PBS -q batch
#PBS -M kh31516@uga.edu 
#PBS -m ae

################################################################################################################################
###### complete pipeline (picard, GATK, annova, Mutect, Mutect2, Germline mutation, Depth of coverage, Callablebases) ##########
##### WGS analysis needs to change the interval region #####

result='/scratch/kh31516/Original_Mammary/results/CMT-2' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Original_Mammary/data/SRR7780922'
dataT='/scratch/kh31516/Original_Mammary/data/SRR7780923'
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/Lab_shared_PanCancer/script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


####### Download #######
module load SRA-Toolkit/2.9.1-centos_linux64

mkdir -p ${result}
mkdir -p $dataN
cd $dataN

fastq-dump --split-files --gzip SRR7780922
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780922 $sraRunTable
mkdir -p $dataT
cd $dataT

fastq-dump --split-files --gzip SRR7780923
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780923 $sraRunTable

####### BWA Mapping #######
cd ${result}
module load BWA/0.7.17-foss-2016b

# Tumor
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR7780923_1.fastq.gz > ${result}/SRR7780923_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR7780923_2.fastq.gz > ${result}/SRR7780923_2.sai

time bwa sampe $reference/canFam3.fa \
${result}/SRR7780923_1.sai \
${result}/SRR7780923_2.sai \
$dataT/SRR7780923_1.fastq.gz \
$dataT/SRR7780923_2.fastq.gz \
> ${result}/SRR7780923.sam

# Normal
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR7780922_1.fastq.gz > ${result}/SRR7780922_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR7780922_2.fastq.gz > ${result}/SRR7780922_2.sai

time bwa sampe $reference/canFam3.fa \
${result}/SRR7780922_1.sai \
${result}/SRR7780922_2.sai \
$dataN/SRR7780922_1.fastq.gz \
$dataN/SRR7780922_2.fastq.gz \
> ${result}/SRR7780922.sam


####### Convert sam to bam file #######

module load SAMtools/1.9-foss-2016b

samtools view -bS ${result}/SRR7780923.sam > ${result}/SRR7780923.bam
samtools view -bS ${result}/SRR7780922.sam > ${result}/SRR7780922.bam

# get header information
#grep '@SQ\|@PG' ${result}/SRR7780922.sam > ${result}/header_N
#grep '@SQ\|@PG' ${result}/SRR7780923.sam > ${result}/header_T
samtools view -H ${source_data}/SRR7780922.bam > ${result}/header_N
samtools view -H ${source_data}/SRR7780923.bam > ${result}/header_T

# exclude unmapped reads based on FLAG
cat ${result}/SRR7780922.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${result}/SRR7780922-cleaned.sam
cat ${result}/header_N ${result}/SRR7780922-cleaned.sam > ${result}/fooN
mv ${result}/fooN ${result}/SRR7780922-cleaned.sam

cat ${result}/SRR7780923.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${result}/SRR7780923-cleaned.sam
cat ${result}/header_T ${result}/SRR7780923-cleaned.sam > ${result}/fooT
mv ${result}/fooT ${result}/SRR7780923-cleaned.sam

####### Picard #######
# picard sort
cd ${result}
module load picard/2.16.0-Java-1.8.0_144

java -Xmx32g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${result}/SRR7780922-cleaned.sam \
O=${result}/SRR7780922_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=SRR7780922_normal

java -Xmx32g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=${result}/SRR7780923-cleaned.sam \
O=${result}/SRR7780923_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=SRR7780923_tumor

# picard MarkDuplicates
java -Xmx32g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${result}/SRR7780922_rg_added_sorted.bam \
O=${result}/SRR7780922_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${result}/SRR7780922-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


java -Xmx32g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=${result}/SRR7780923_rg_added_sorted.bam \
O=${result}/SRR7780923_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${result}/SRR7780923-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### GATK Realign #######
# Generating interval file for sort.bam
module load GATK/3.8-1-Java-1.8.0_144
java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa \
-I ${result}/SRR7780922_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${result}/SRR7780922_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa \
-I ${result}/SRR7780923_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${result}/SRR7780923_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa \
-I ${result}/SRR7780922_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${result}/SRR7780922_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa \
-I ${result}/SRR7780923_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${result}/SRR7780923_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 

cd ${result}
module load MuTect/1.1.7-Java-1.7.0_80

time java -Xmx32g -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals $reference/Canis_familiaris.CanFam3.interval_list \
--input_file:normal ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${result}/CMT-2_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${result}/CMT-2_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf > ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${result}/CMT-2_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${result}/CMT-2_PASS.stat

python $script/Filter_MutectStat_5steps.py ${result}/CMT-2_PASS.stat ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


####### Mutect2 #######

cd ${result}
module load SAMtools/1.9-foss-2016b
module load GATK/4.1.0.0-Java-1.8.0_144

Normal_sample=$(samtools view -H ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R $reference/canFam3.fa \
-I ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
-L $reference/Canis_familiaris.CanFam3.interval_list \
-O ${result}/CMT-2_MuTect2_GATK4_noDBSNP.vcf


###### Germline mutation preparation ######
################ GATK 3 ###################

module load GATK/3.8-1-Java-1.8.0_144
cd ${result}

# Variant calling
java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I \
${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I \
${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa \
-V ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa \
-V ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation
module load Perl/5.26.1-GCCcore-6.4.0

perl $annovar_index/convert2annovar.pl -format vcf4old ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output
python $script/Add_GeneName_N_Signature.py ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python $script/Add_GeneName_N_Signature.py ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


#### GATK DepthofCoverage ####
cd ${result}
module load GATK/3.8-1-Java-1.8.0_144

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${result}/SRR7780922_DepthofCoverage_CDS.bed 

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${result}/SRR7780923_DepthofCoverage_CDS.bed 


#### GATK callable ####
cd ${result}
module load GATK/3.8-1-Java-1.8.0_144

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
 -R $reference/canFam3.fa \
 -I SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
 -summary SRR7780922_table.txt \
 -o SRR7780922_callable_status.bed 

java -Xmx32g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
 -R $reference/canFam3.fa \
 -I SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
 -summary SRR7780923_table.txt \
 -o SRR7780923_callable_status.bed 

####### remove unneeded files #######
rm ${result}/*_sorted_dedupped_removed.bam
rm ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS ${result}/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm ${result}/*.bai
rm ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${result}/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${result}/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${result}/SRR7780922.sam ${result}/SRR7780923.sam 
rm ${result}/header_N ${result}/header_T
rm ${result}/SRR7780922_rg_added_sorted.bam ${result}/SRR7780923_rg_added_sorted.bam
rm ${result}/SRR7780922-cleaned.sam ${result}/SRR7780923-cleaned.sam
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm ${result}/*.sai

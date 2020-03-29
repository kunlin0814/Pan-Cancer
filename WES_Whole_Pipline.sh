#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=360:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -q batch
#PBS -M kh31516@uga.edu 
#PBS -m ae
###### complete pipeline (including callable base) ##########


result='/scratch/kh31516/Original_Mammary/results/CMT-2' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
dataN='/scratch/kh31516/Original_Mammary/data/SRR7780922'
dataT='/scratch/kh31516/Original_Mammary/data/SRR7780923'
annovar_index='/scratch/kh31516/Melanoma/Melanoma_source/annovar'
script='/scratch/kh31516/variant_calling/scripts'
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0


####### Download #######
mkdir $result
mkdir $dataN
cd $dataN

fastq-dump --split-files SRR7780922
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780922 $sraRunTable
mkdir $dataT
cd $dataT

fastq-dump --split-files SRR7780923
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780923 $sraRunTable


####### BWA Mapping #######
cd $result
# Tumor
time bwa aln $reference/canFam3.fa $dataT/SRR7780923_1.fastq > $result/SRR7780923_1.sai
time bwa aln $reference/canFam3.fa $dataT/SRR7780923_2.fastq > $result/SRR7780923_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780923_1.sai \
$result/SRR7780923_2.sai \
$dataT/SRR7780923_1.fastq \
$dataT/SRR7780923_2.fastq \
> $result/SRR7780923.sam

# Normal
time bwa aln $reference/canFam3.fa $dataN/SRR7780922_1.fastq > $result/SRR7780922_1.sai
time bwa aln $reference/canFam3.fa $dataN/SRR7780922_2.fastq > $result/SRR7780922_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780922_1.sai \
$result/SRR7780922_2.sai \
$dataN/SRR7780922_1.fastq \
$dataN/SRR7780922_2.fastq \
> $result/SRR7780922.sam




####### Convert sam to bam file #######
samtools view -bS $result/SRR7780923.sam > $result/SRR7780923.bam
samtools view -bS $result/SRR7780922.sam > $result/SRR7780922.bam




####### Delete #######
rm $dataN/*.fastq $dataT/*.fastq
rm $result/*.sai




####### Picard #######
# get header information
grep '@SQ\|@PG' $result/SRR7780922.sam > $result/header_N
grep '@SQ\|@PG' $result/SRR7780923.sam > $result/header_T

# exclude unmapped reads based on FLAG
cat $result/SRR7780922.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR7780922-cleaned.sam
cat $result/header_N $result/SRR7780922-cleaned.sam > $result/fooN
mv $result/fooN $result/SRR7780922-cleaned.sam
cat $result/SRR7780923.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR7780923-cleaned.sam
cat $result/header_T $result/SRR7780923-cleaned.sam > $result/fooT
mv $result/fooT $result/SRR7780923-cleaned.sam

# picard sort
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR7780922-cleaned.sam O=$result/SRR7780922_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR7780923-cleaned.sam O=$result/SRR7780923_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR7780922_rg_added_sorted.bam O=$result/SRR7780922_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SRR7780922-output.metrics REMOVE_SEQUENCING_DUPLICATES=true
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR7780923_rg_added_sorted.bam O=$result/SRR7780923_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SRR7780923-output.metrics REMOVE_SEQUENCING_DUPLICATES=true





####### Remove unneeded files #######
rm $result/SRR7780922.sam $result/SRR7780923.sam 
rm $result/header_N $result/header_T
rm $result/SRR7780922_rg_added_sorted.bam $result/SRR7780923_rg_added_sorted.bam
rm $result/SRR7780922-cleaned.sam $result/SRR7780923-cleaned.sam



####### GATK Realign #######
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SRR7780922_rg_added_sorted_dedupped_removed.bam -o $result/SRR7780922_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SRR7780923_rg_added_sorted_dedupped_removed.bam -o $result/SRR7780923_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

# realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SRR7780922_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SRR7780922_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SRR7780923_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SRR7780923_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores



######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/CMT-2_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/CMT-2_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf





####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference


####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
rm $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/CMT-2_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm $result/*.bai

###### Germline mutation preparation ######

cd $result
module load SAMtools/1.9-foss-2016b
module load GATK/3.8-1-Java-1.8.0_144

################ GATK Annovar ############
################ GATK ###################
cd $result
samtools index $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam
samtools index $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam

# Variant calling
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.vcf
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa -V $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa -V $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS
awk '$7 == "PASS" {print $0}' $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf > $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS
# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS > $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
perl $annovar_index/convert2annovar.pl -format vcf4old $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS > $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index
# add gene names to annovar output
python $script/Add_GeneName_N_Signature.py $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.96.gtf_geneNamePair.txt
python $script/Add_GeneName_N_Signature.py $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.96.gtf_geneNamePair.txt


############### Delete files ##############
rm $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput


module load GATK/3.8-1-Java-1.8.0_144

cd $result


#### GATK DepthofCoverage ####

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o $result/SRR7780922_DepthofCoverage_CDS.bed 

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o $result/SRR7780923_DepthofCoverage_CDS.bed 


#### GATK callable ####

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
 -R $reference/canFam3.fa \
 -I SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
 -summary SRR7780976_table.txt \
 -o SRR7780922_callable_status.bed 

time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
 -R $reference/canFam3.fa \
 -I SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
 -summary SRR7780976_table.txt \
 -o SRR7780923_callable_status.bed 

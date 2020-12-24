#!/bin/bash
#PBS -N Trimmomatic-i_03A6-WES-20
#PBS -l walltime=170:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=60gb
#PBS -q batch
###### complete pipeline (including callable base) ##########


result='/scratch/kh31516/Pan_cancer/glioma/results/Trimmomatic_WES/BQ20/results/i_03A6' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Pan_cancer/glioma/results/Trimmomatic_WES/data/SRR10351810'
dataT='/scratch/kh31516/Pan_cancer/glioma/results/Trimmomatic_WES/data/SRR10351811'
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/Lab_shared_PanCancer/script'
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0


####### Download #######
#mkdir $result
#mkdir $dataN
#cd $dataN

#fastq-dump --gzip --split-files SRR10351810
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR10351810 $sraRunTable
#mkdir $dataT
#cd $dataT

#fastq-dump --gzip --split-files SRR10351811
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR10351811 $sraRunTable


mkdir $result
cd $result

###  Trimmomatic

module load Trimmomatic/0.36-Java-1.8.0_144

time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 \
$dataN/SRR10351810_1.fastq.gz $dataN/SRR10351810_2.fastq.gz \
$dataN/SRR10351810_lane1_forward_paired.fq.gz $dataN/SRR10351810_lane1_forward_unpaired.fq.gz \
$dataN/SRR10 351810_lane1_reverse_paired.fq.gz $dataN/SRR10351810_lane1_reverse_unpaired.fq.gz \
ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36



time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 \
$dataT/SRR10351811_1.fastq.gz $dataT/SRR10351811_2.fastq.gz \
$dataT/SRR10351811_lane1_forward_paired.fq.gz $dataT/SRR10351811_lane1_forward_unpaired.fq.gz \
$dataT/SRR10351811_lane1_reverse_paired.fq.gz $dataT/SRR10351811_lane1_reverse_unpaired.fq.gz \
ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

####### BWA Mapping #######
cd $result
# Tumor
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR10351811_lane1_forward_paired.fq.gz > $result/SRR10351811_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR10351811_lane1_reverse_paired.fq.gz > $result/SRR10351811_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR10351811_1.sai \
$result/SRR10351811_2.sai \
$dataT/SRR10351811_lane1_forward_paired.fq.gz \
$dataT/SRR10351811_lane1_reverse_paired.fq.gz \
> $result/SRR10351811.sam

# Normal
time bwa aln -t 4 $reference/canFam3.fa $dataN//SRR10351810_lane1_forward_paired.fq.gz > $result/SRR10351810_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR10351810_lane1_reverse_paired.fq.gz > $result/SRR10351810_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR10351810_1.sai \
$result/SRR10351810_2.sai \
$dataN//SRR10351810_lane1_forward_paired.fq.gz \
$dataN/SRR10351810_lane1_reverse_paired.fq.gz \
> $result/SRR10351810.sam




####### Convert sam to bam file #######
samtools view -bS $result/SRR10351811.sam > $result/SRR10351811.bam
samtools view -bS $result/SRR10351810.sam > $result/SRR10351810.bam




####### Delete #######
#rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm $result/*.sai


####### Picard #######
# get header information
cd $result
grep '@SQ\|@PG' $result/SRR10351810.sam > $result/header_N
grep '@SQ\|@PG' $result/SRR10351811.sam > $result/header_T

# exclude unmapped reads based on FLAG
cat $result/SRR10351810.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR10351810-cleaned.sam
cat $result/header_N $result/SRR10351810-cleaned.sam > $result/fooN
mv $result/fooN $result/SRR10351810-cleaned.sam
cat $result/SRR10351811.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR10351811-cleaned.sam
cat $result/header_T $result/SRR10351811-cleaned.sam > $result/fooT
mv $result/fooT $result/SRR10351811-cleaned.sam

# picard sort
java -Xmx16g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR10351810-cleaned.sam \
O=$result/SRR10351810_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

java -Xmx16g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR10351811-cleaned.sam \
O=$result/SRR10351811_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
java -Xmx16g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR10351810_rg_added_sorted.bam \
O=$result/SRR10351810_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=$result/SRR10351810-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


java -Xmx16g -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR10351811_rg_added_sorted.bam \
O=$result/SRR10351811_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=$result/SRR10351811-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### Remove unneeded files #######
rm $result/SRR10351810.sam $result/SRR10351811.sam
rm $result/header_N $result/header_T
rm $result/SRR10351810_rg_added_sorted.bam $result/SRR10351811_rg_added_sorted.bam
rm $result/SRR10351810-cleaned.sam $result/SRR10351811-cleaned.sam


####### GATK Realign #######
# Generating interval file for sort.bam
cd $result
java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa \
-I $result/SRR10351810_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o $result/SRR10351810_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa \
-I $result/SRR10351811_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o $result/SRR10351811_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa \
-I $result/SRR10351810_rg_added_sorted_dedupped_removed.bam \
-targetIntervals $result/SRR10351810_rg_added_sorted_dedupped_removed.bam.intervals \
-o $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa \
-I $result/SRR10351811_rg_added_sorted_dedupped_removed.bam \
-targetIntervals $result/SRR10351811_rg_added_sorted_dedupped_removed.bam.intervals \
-o $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores


######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
cd $result
module load MuTect/1.1.7-Java-1.7.0_80

java -Xmx16g -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbSNP_2020/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/i_03A6_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/i_03A6_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf


####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
> $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# add gene name
python $script/Add_GeneName_N_Signature.py $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


############ 5 Steps filtering ###################
grep -w KEEP $result/i_03A6_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > $result/i_03A6_PASS.stat
python $script/Filter_MutectStat_5steps.py $result/i_03A6_PASS.stat $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $reference/annovar_CanFam3.1.99.gtf

# add gene name
python $script/Add_GeneName_N_Signature.py $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

##################################################



####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
rm $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm $result/*.bai


####### Mutect2 #######
#java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/SRR10351618_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf



###### Germline mutation preparation ######

cd $result
module load SAMtools/1.9-foss-2016b
module load GATK/3.8-1-Java-1.8.0_144

################ GATK Annovar ############
################ GATK ###################

cd $result
samtools index $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam
samtools index $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam

# Variant calling
java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I \
$result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference/canFam3.fa -I \
$result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa \
-V $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $reference/canFam3.fa \
-V $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation

perl $annovar_index/convert2annovar.pl -format vcf4old $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

#annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output
python $script/Add_GeneName_N_Signature.py $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python $script/Add_GeneName_N_Signature.py $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

############### Delete files ##############
rm $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

module load GATK/3.8-1-Java-1.8.0_144

cd $result

#### GATK DepthofCoverage ####

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I $result/SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o $result/SRR10351810_DepthofCoverage_CDS.bed 

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R $reference/canFam3.fa \
-I $result/SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o $result/SRR10351811_DepthofCoverage_CDS.bed 


#### GATK callable ####

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R $reference/canFam3.fa \
-I SRR10351810_rg_added_sorted_dedupped_removed.realigned.bam \
-summary SRR10351810_table.txt \
-o SRR10351810_callable_status.bed 

java -Xmx16g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R $reference/canFam3.fa \
-I SRR10351811_rg_added_sorted_dedupped_removed.realigned.bam \
-summary SRR10351811_table.txt \
-o SRR10351811_callable_status.bed
  
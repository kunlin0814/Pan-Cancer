#!/bin/bash
#PBS -N CMT-263
#PBS -l walltime=128:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=60gb
#PBS -q batch

result='/scratch/yf94402/Canine_Mammary_cancer/results/WXS/CMT-263'
reference='/scratch/yf94402/variant_calling/source'
dataN='/scratch/yf94402/Canine_Mammary_cancer/WXS/SRR7780735'
dataT='/scratch/yf94402/Canine_Mammary_cancer/WXS/SRR7780736'

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
fastq-dump --split-files --gzip SRR7780735
mkdir $dataT
cd $dataT
fastq-dump --split-files --gzip SRR7780736




####### BWA Mapping #######
cd $result
# Tumor
time bwa aln $reference/canFam3.fa $dataT/SRR7780736_1.fastq.gz > $result/SRR7780736_1.sai
time bwa aln $reference/canFam3.fa $dataT/SRR7780736_2.fastq.gz > $result/SRR7780736_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780736_1.sai \
$result/SRR7780736_2.sai \
$dataT/SRR7780736_1.fastq.gz \
$dataT/SRR7780736_2.fastq.gz \
> $result/SRR7780736.sam

# Normal
time bwa aln $reference/canFam3.fa $dataN/SRR7780735_1.fastq.gz > $result/SRR7780735_1.sai
time bwa aln $reference/canFam3.fa $dataN/SRR7780735_2.fastq.gz > $result/SRR7780735_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780735_1.sai \
$result/SRR7780735_2.sai \
$dataN/SRR7780735_1.fastq.gz \
$dataN/SRR7780735_2.fastq.gz \
> $result/SRR7780735.sam




####### Convert sam to bam file #######
samtools view -bS $result/SRR7780736.sam > $result/SRR7780736.bam
samtools view -bS $result/SRR7780735.sam > $result/SRR7780735.bam




####### Delete #######
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm $result/*.sai




####### Picard #######
# get header information
grep '@SQ\|@PG' $result/SRR7780735.sam > $result/header_N
grep '@SQ\|@PG' $result/SRR7780736.sam > $result/header_T

# exclude unmapped reads based on FLAG
cat $result/SRR7780735.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR7780735-cleaned.sam
cat $result/header_N $result/SRR7780735-cleaned.sam > $result/fooN
mv $result/fooN $result/SRR7780735-cleaned.sam
cat $result/SRR7780736.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SRR7780736-cleaned.sam
cat $result/header_T $result/SRR7780736-cleaned.sam > $result/fooT
mv $result/fooT $result/SRR7780736-cleaned.sam

# picard sort
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR7780735-cleaned.sam O=$result/SRR7780735_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SRR7780736-cleaned.sam O=$result/SRR7780736_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR7780735_rg_added_sorted.bam O=$result/SRR7780735_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SRR7780735-output.metrics REMOVE_SEQUENCING_DUPLICATES=true
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SRR7780736_rg_added_sorted.bam O=$result/SRR7780736_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SRR7780736-output.metrics REMOVE_SEQUENCING_DUPLICATES=true





####### Remove unneeded files #######
rm $result/SRR7780735.sam $result/SRR7780736.sam 
rm $result/header_N $result/header_T
rm $result/SRR7780735_rg_added_sorted.bam $result/SRR7780736_rg_added_sorted.bam
rm $result/SRR7780735-cleaned.sam $result/SRR7780736-cleaned.sam





####### GATK Realign #######
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SRR7780735_rg_added_sorted_dedupped_removed.bam -o $result/SRR7780735_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SRR7780736_rg_added_sorted_dedupped_removed.bam -o $result/SRR7780736_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

# realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SRR7780735_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SRR7780735_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SRR7780735_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SRR7780736_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SRR7780736_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SRR7780736_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores





######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal $result/SRR7780735_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/SRR7780736_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/CMT-263_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/CMT-263_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/SRR7780736_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/SRR7780735_rg_added_sorted_dedupped_removed.realigned.bam \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf





####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference






####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
rm $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm $result/*.bai



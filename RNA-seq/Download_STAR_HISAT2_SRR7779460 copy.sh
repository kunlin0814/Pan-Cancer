#!/bin/bash
#PBS -N SRR7779460
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=128:00:00
#PBS -l mem=60gb

### INDEX the genome for STAR and HISAT2 respectively first ###

genomeDir='/scratch/kh31516/Melanoma/RNA-seq/STAR'
data='/scratch/kh31516/Melanoma/RNA-seq/data/SRR7779460'
result='/scratch/kh31516/Melanoma/RNA-seq/results/SRR7779460'
hisat_source='/scratch/kh31516/Melanoma/RNA-seq/HISAT2'
source='/scratch/kh31516/Melanoma_source'

mkdir $result
mkdir $result/STAR
mkdir $result/HISAT2

ml STAR/2.6.1c-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load SAMtools/1.9-foss-2016b
module load HISAT2/2.1.0-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
module load Cufflinks/2.2.1-foss-2016b

####### Download #######
mkdir $data
cd $data
fastq-dump --gzip SRR7779460






############### Germline Muation Detection with STAR ##################
####### STAR #######
# 1pass
cd $result/STAR
mkdir 1pass
cd $result/STAR/1pass
STAR --genomeDir $genomeDir --readFilesIn $data/SRR7779460.fastq.gz --readFilesCommand zcat --runThreadN 4
# reindex reference genome with junction info from 1pass
genomeDir2=$result/STAR/canFam3_2pass
mkdir $genomeDir2
cp $source/canFam3.fa $genomeDir2
cd $genomeDir2
STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles canFam3.fa --sjdbFileChrStartEnd $result/STAR/1pass/SJ.out.tab --sjdbOverhang 75 --runThreadN 24
# 2pass with new reference index
cd $result/STAR
mkdir 2pass
cd 2pass
STAR --genomeDir $genomeDir2 --readFilesIn $data/SRR7779460.fastq.gz --readFilesCommand zcat --runThreadN 4

cd $result/STAR
samtools view -bS $result/STAR/2pass/Aligned.out.sam > $result/STAR/SRR7779460.bam
rm -r $result/STAR/1pass $result/STAR/canFam3_2pass $result/STAR/2pass

####### Picard #######
# Sort the bam file
module load picard/2.16.0-Java-1.8.0_144
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/STAR/SRR7779460.bam O=$result/STAR/SRR7779460-rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
# Mark duplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/STAR/SRR7779460-rg_added_sorted.bam O=$result/STAR/SRR7779460_rg_added_sorted_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=SRR7779460-output.metrics

####### GATK #######
#Split'N'Trim and reassign mapping qualities
module load GATK/3.8-1-Java-1.8.0_144
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R $source/canFam3.fa -I $result/STAR/SRR7779460_rg_added_sorted_dedupped.bam -o $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $source/canFam3.fa -I $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam -o $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam.intervals --allow_potentially_misencoded_quality_scores
# Realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $source/canFam3.fa -I $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam -targetIntervals $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam.intervals -o $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam --allow_potentially_misencoded_quality_scores
# Variant calling
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R $source/canFam3.fa -I $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam.vcf
# Variant filtering
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R $source/canFam3.fa -V $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf






############### Expression Level with HISAT2 #########################
####### HISAT2 #######
module load HISAT2/2.1.0-foss-2016b
cd $result/HISAT2
hisat2 -p 8 -x $hisat_source/canFam3 \
-U $data/SRR7779460.fastq.gz \
-S $result/HISAT2/SRR7779460.sam

####### sort bam file #######
module load SAMtools/1.9-foss-2016b
samtools view -bS $result/HISAT2/SRR7779460.sam > $result/HISAT2/SRR7779460.bam
mkdir $result/HISAT2/tmp
module load SAMtools/1.9-foss-2016b
samtools sort -T $result/HISAT2/tmp/ -o $result/HISAT2/SRR7779460-sorted.bam $result/HISAT2/SRR7779460.bam
rm -r $result/HISAT2/tmp
rm $result/HISAT2/SRR7779460.sam

####### cufflinks #######
cufflinks -p 8 -g $source/canFam3.gtf -b $source/canFam3.fa -u -o $result/HISAT2/ $result/HISAT2/SRR7779460-sorted.bam






####### Delete #######
rm $data/*.fastq.gz
rm $result/STAR/SRR7779460-rg_added_sorted.bam $result/STAR/SRR7779460_rg_added_sorted_dedupped.bam $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam
rm $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.bam.intervals $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam
rm $result/STAR/*.bai $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam.vcf $result/STAR/SRR7779460-rg_added_sorted_dedupped_split.realigned.bam.vcf.idx
rm $result/HISAT2/SRR7779460-sorted.bam

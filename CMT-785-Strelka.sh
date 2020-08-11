#!/bin/bash
#PBS -N CMT-785
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -q batch

###### complete pipeline (including callable) ##########


##data='/scratch/kh31516/Original_Mammary/data'
result='/scratch/kh31516/Original_Mammary/results/CMT-785' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
dataN='/scratch/kh31516/Original_Mammary/data/SRR7780807'
dataT='/scratch/kh31516/Original_Mammary/data/SRR7780755'
annovar_index='/scratch/kh31516/Melanoma/Melanoma_source/annovar'
script='/scratch/kh31516/variant_calling/scripts'
strelka_result='/scratch/kh31516/Original_Mammary/strelka_result'
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0

: "
####### Download #######
mkdir $result
mkdir $dataN
cd $dataN
#fasterq-dump SRR7780807 -t /scratch/kh31516/Gtex/tmp/ -e 4 --split-files
fastq-dump --split-files SRR7780807
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780807 $sraRunTable
mkdir $dataT
cd $dataT
#fasterq-dump SRR7780755 -t /scratch/kh31516/Gtex/tmp/ -e 4 --split-files
fastq-dump --split-files SRR7780755
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh SRR7780755 $sraRunTable


####### BWA Mapping #######
cd $result
# Tumor
time bwa aln $reference/canFam3.fa $dataT/SRR7780755_1.fastq > $result/SRR7780755_1.sai
time bwa aln $reference/canFam3.fa $dataT/SRR7780755_2.fastq > $result/SRR7780755_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780755_1.sai \
$result/SRR7780755_2.sai \
$dataT/SRR7780755_1.fastq \
$dataT/SRR7780755_2.fastq \
> $result/SRR7780755.sam

# Normal
time bwa aln $reference/canFam3.fa $dataN/SRR7780807_1.fastq > $result/SRR7780807_1.sai
time bwa aln $reference/canFam3.fa $dataN/SRR7780807_2.fastq > $result/SRR7780807_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR7780807_1.sai \
$result/SRR7780807_2.sai \
$dataN/SRR7780807_1.fastq \
$dataN/SRR7780807_2.fastq \
> $result/SRR7780807.sam




####### Convert sam to bam file #######
samtools view -bS $result/SRR7780755.sam > $result/SRR7780755.bam
samtools view -bS $result/SRR7780807.sam > $result/SRR7780807.bam
"
strelka_folder='/scratch/kh31516/Pan_cancer'
cd $result

samtools index $result/SRR7780807_rg_added_sorted_dedupped_removed.realigned.bam
samtools index $result/SRR7780755_rg_added_sorted_dedupped_removed.realigned.bam

$strelka_folder/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $result/SRR7780807_rg_added_sorted_dedupped_removed.realigned.bam \
    --tumorBam $result/SRR7780755_rg_added_sorted_dedupped_removed.realigned.bam \
    --referenceFasta $reference/canFam3.fa \
    --runDir $strelka_result/CMT-785-demo_somatic

# your result will be generated in demo_somatic/results/variants
$strelka_result/CMT-785-demo_somatic/runWorkflow.py -m local -j 20
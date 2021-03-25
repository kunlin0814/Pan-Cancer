#!/bin/bash
#SBATCH --job-name=re-run_DD0001a_WES         # Job name (DD0001_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=8           # Number of CPU cores per task (4)
#SBATCH --mem=40G                   # Job memory limit (10 GB)
#SBATCH --time=50:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=DD0001a_WES.%j.out    # Standard output log
#SBATCH --error=DD0001a_WES.%j.err     # Standard error log

thread=8
Bioproject='PRJNA552034'
Normal_Run='ERR1681499'
Tumor_Run='ERR1681500'
SampleName='DD0001a'
results='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/results/'${SampleName} #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/data/'${SampleName}'/'${Normal_Run}
dataT='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/data/'${SampleName}'/'${Tumor_Run}
MuTect_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/Mutect/'${SampleName}''
# MuTect2_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/Mutect2/'${SampleName}''
# MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
# Germline_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/Mutect/Germline/'${SampleName}''
# DepthOfCoverage='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/Mutect/DepthOfCoverage/'${SampleName}''
strelka_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_whole/Mutect/strelka/'${SampleName}''
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/home/kh31516/kh31516_Lab_Share_script'



####### Download #######
module load SRA-Toolkit/2.9.6-1-centos_linux64
ml Anaconda3/2020.02
mkdir -p $dataT
mkdir -p $dataN
mkdir -p ${results}
mkdir -p ${MuTect_out}
# mkdir -p ${MuTect2_out}
# mkdir -p ${Germline_out}
# mkdir -p ${DepthOfCoverage}
# mkdir -p $strelka_out


cd $dataN

fastq-dump --split-files --gzip ${Normal_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Normal_Run} $sraRunTable

cd $dataT

fastq-dump --split-files --gzip ${Tumor_Run}
#/scratch/kh31516/Pan_cancer/PRJNA247493/download_sra_sample.sh ${Tumor_Run} $sraRunTable

####### BWA Mapping #######
cd ${results}
ml BWA/0.7.17-GCC-8.3.0

# Tumor
time bwa aln -t ${thread} ${reference}/canFam3.fa $dataT/${Tumor_Run}_1.fastq.gz > ${results}/${Tumor_Run}_1.sai
time bwa aln -t ${thread} ${reference}/canFam3.fa $dataT/${Tumor_Run}_2.fastq.gz > ${results}/${Tumor_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Tumor_Run}_1.sai \
${results}/${Tumor_Run}_2.sai \
$dataT/${Tumor_Run}_1.fastq.gz \
$dataT/${Tumor_Run}_2.fastq.gz \
> ${results}/${Tumor_Run}.sam

# Normal
time bwa aln -t ${thread} ${reference}/canFam3.fa $dataN/${Normal_Run}_1.fastq.gz > ${results}/${Normal_Run}_1.sai
time bwa aln -t ${thread} ${reference}/canFam3.fa $dataN/${Normal_Run}_2.fastq.gz > ${results}/${Normal_Run}_2.sai

time bwa sampe ${reference}/canFam3.fa \
${results}/${Normal_Run}_1.sai \
${results}/${Normal_Run}_2.sai \
$dataN/${Normal_Run}_1.fastq.gz \
$dataN/${Normal_Run}_2.fastq.gz \
> ${results}/${Normal_Run}.sam


# ####### Convert sam to bam file #######

# cd ${results}
# module load SAMtools/1.9-GCC-8.3.0

# samtools view -bS ${results}/${Tumor_Run}.sam > ${results}/${Tumor_Run}.bam
# samtools view -bS ${results}/${Normal_Run}.sam > ${results}/${Normal_Run}.bam

# # get header information
# #grep '@SQ\|@PG' ${results}/${Normal_Run}.sam > ${results}/header_N
# #grep '@SQ\|@PG' ${results}/${Tumor_Run}.sam > ${results}/header_T
# samtools view -H ${results}/${Normal_Run}.bam > ${results}/header_N
# samtools view -H ${results}/${Tumor_Run}.bam > ${results}/header_T

# # exclude unmapped reads based on FLAG
# cat ${results}/${Normal_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Normal_Run}-cleaned.sam
# cat ${results}/header_N ${results}/${Normal_Run}-cleaned.sam > ${results}/fooN
# mv ${results}/fooN ${results}/${Normal_Run}-cleaned.sam

# cat ${results}/${Tumor_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Tumor_Run}-cleaned.sam
# cat ${results}/header_T ${results}/${Tumor_Run}-cleaned.sam > ${results}/fooT
# mv ${results}/fooT ${results}/${Tumor_Run}-cleaned.sam
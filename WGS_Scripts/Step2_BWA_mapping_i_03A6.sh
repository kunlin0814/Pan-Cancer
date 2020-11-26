#!/bin/bash
#SBATCH --job-name=i_03A6         # Job name (i_03A6)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=8           # Number of CPU co4es per task (4)
#SBATCH --mem=70G                   # Job memory limit (10 GB)
#SBATCH --time=166:00:00              # Time 1it hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=i_03A6.%j.out    # Standard output log
#SBATCH --error=i_03A6.%j.err     # Standard1ror4log


sample_name='i_03A6'
Normal_Run='SRR10362478'
Tumor_Run='SRR10362479'
results='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/results/'${sample_name} #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/'${sample_name}/${Normal_Run}
dataT='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/'${sample_name}/${Tumor_Run}
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/home/kh31516/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect/'${sample_name}
MuTect2_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect2/'${sample_name}
Germline_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Germline/'${sample_name}
DepthOfCoverage='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/DepthOfCoverage/'${sample_name}
strelka_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/strelka/'

####### BWA Mapping #######
cd ${results}
ml BWA/0.7.17-GCC-8.3.0

# Tumor
bwa aln -t 8 ${reference}/canFam3.fa $dataT/${Tumor_Run}_1.fastq.gz > ${results}/${Tumor_Run}_1.sai
bwa aln -t 8 ${reference}/canFam3.fa $dataT/${Tumor_Run}_2.fastq.gz > ${results}/${Tumor_Run}_2.sai

bwa sampe ${reference}/canFam3.fa \
${results}/${Tumor_Run}_1.sai \
${results}/${Tumor_Run}_2.sai \
$dataT/${Tumor_Run}_1.fastq.gz \
$dataT/${Tumor_Run}_2.fastq.gz \
> ${results}/${Tumor_Run}.sam

# Normal
bwa aln -t 8 ${reference}/canFam3.fa $dataN/${Normal_Run}_1.fastq.gz > ${results}/${Normal_Run}_1.sai
bwa aln -t 8 ${reference}/canFam3.fa $dataN/${Normal_Run}_2.fastq.gz > ${results}/${Normal_Run}_2.sai

bwa sampe ${reference}/canFam3.fa \
${results}/${Normal_Run}_1.sai \
${results}/${Normal_Run}_2.sai \
$dataN/${Normal_Run}_1.fastq.gz \
$dataN/${Normal_Run}_2.fastq.gz \
> ${results}/${Normal_Run}.sam


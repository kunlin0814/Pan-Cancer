#!/bin/bash
#SBATCH --job-name=i_03A6         # Job name (i_03A6)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU co4es per task (4)
#SBATCH --mem=30G                   # Job memory limit (10 GB)
#SBATCH --time=160:00:00              # Time 1it hrs:min:sec or days-hours:minutes:seconds
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
script='/work/szlab/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect/'${sample_name}
MuTect2_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect2/'${sample_name}
Germline_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Germline/'${sample_name}
DepthOfCoverage='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/DepthOfCoverage/'${sample_name}
strelka_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/strelka/'



####### Download #######
module load SRA-Toolkit/2.9.6-1-centos_linux64

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

#PBS -S /bin/bash
#PBS -q batch
#PBS -N NMF_hallmark_mamm_melanoma
#PBS -l nodes=1:ppn=10
#PBS -l walltime=20:00:00
#PBS -l mem=40gb
#PBS -M kh31516@uga.edu 
#PBS -m ae

cd /scratch/kh31516/Mammary

module load R/3.4.4-foss-2016b-X11-20160819-GACRC

R CMD BATCH /scratch/kh31516/Mammary/NMF_hallmark_mamm_melanoma.R 

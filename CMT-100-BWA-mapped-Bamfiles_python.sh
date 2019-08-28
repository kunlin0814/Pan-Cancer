#PBS -S /bin/bash
#PBS -q batch
#PBS -N CMT-100
#PBS -l nodes=1:ppn=4
#PBS -l walltime=100:00:00
#PBS -l mem=150gb
#PBS -M kh31516@uga.edu 
#PBS -m ae
cd /scratch/kh31516/Mammary/results
module load SAMtools/0.1.19-foss-2016b
#samtools view /scratch/kh31516/Mammary/results/CMT-100/SRR7780976.bam > /scratch/kh31516/Mammary/results/CMT-100/SRR7780976.sam
samtools view /scratch/kh31516/Mammary/results/CMT-100/SRR7780979.bam > /scratch/kh31516/Mammary/results/CMT-100/SRR7780979.sam

#python /scratch/kh31516/Mammary/Summarize_BWA_sam.py /scratch/kh31516/Mammary/results/CMT-100/SRR7780976.sam SRR7780976
python /scratch/kh31516/Mammary/Summarize_BWA_sam.py /scratch/kh31516/Mammary/results/CMT-100/SRR7780979.sam SRR7780979

#rm /scratch/kh31516/Mammary/results/CMT-100/SRR7780976.sam
rm /scratch/kh31516/Mammary/results/CMT-100/SRR7780979.sam
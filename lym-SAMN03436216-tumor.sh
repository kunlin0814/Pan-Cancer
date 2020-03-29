#PBS -S /bin/bash
#PBS -N lym-SAMN03436216
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=15gb
#PBS -q batch



result='/scratch/kh31516/Pan_cancer/lym/results/tumor/SAMN03436216' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
#data='/scratch/kh31516/Trio-family/Original_Trio/data/DD0001/ERR1681499-ERR1681454'
script='/work/szlab/kh31516_Lab_Share_script'

module load SAMtools/1.9-foss-2016b


## Analyze the mapping result

cd $result/

module load SAMtools/0.1.19-foss-2016b

samtools view SAMN03436216.bam > SAMN03436216.sam

java -cp $script/ Summarized_BWA \
$result/SAMN03436216.sam \
SAMN03436216 \
SAMN03436216_BWA_summarize
#python3 /scratch/kh31516/Trio-family/Original_Trio/Summarize_BWA_sam.py $result/SAMN03436216.sam SAMN03436216

cat $result/SAMN03436216_BWA_summarize >> /scratch/kh31516/Pan_cancer/lym/results/Total_Tumor_Lymphoma_BWA_summarization.txt

rm $result/SAMN03436216.sam
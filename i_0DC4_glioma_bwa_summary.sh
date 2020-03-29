#PBS -S /bin/bash
#PBS -N i_0DC4_glioma
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -q batch



result='/scratch/kh31516/Pan_cancer/glioma/results/i_0DC4' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
#data='/scratch/kh31516/Trio-family/Original_Trio/data/DD0001/ERR1681499-ERR1681454'
script='/work/szlab/kh31516_Lab_Share_script'

module load SAMtools/1.9-foss-2016b


## Analyze the mapping result

cd $result/

samtools view NA.bam > NA.sam # normal
samtools view SRR10351582.bam > SRR10351582.sam # tumor

java -cp $script/ Summarized_BWA \
$result/NA.sam \
NA \
NA_BWA_summarize

java -cp $script/ Summarized_BWA \
$result/SRR10351582.sam \
SRR10351582 \
SRR10351582_BWA_summarize
#python3 /scratch/kh31516/Trio-family/Original_Trio/Summarize_BWA_sam.py $result/SAMN03436273.sam SAMN03436273
cat $result/NA_BWA_summarize >> /scratch/kh31516/Pan_cancer/glioma/Total_Normal_glioma_BWA_summarization.txt
cat $result/SRR10351582_BWA_summarize >> /scratch/kh31516/Pan_cancer/glioma/Total_Tumor_glioma_BWA_summarization.txt

rm $result/SAMN03436273.sam
rm $result/NA.sam
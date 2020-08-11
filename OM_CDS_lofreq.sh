#!/bin/bash
#PBS -q batch
#PBS -N OM-CDS-lofreq
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -l mem=60gb


reference='/work/szlab/Lab_shared_PanCancer/source'
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
Lofreq_folder='/scratch/kh31516/Original_Melanoma/store/lofreq'
VarScan2_folder='/scratch/kh31516/Pan_cancer/lymphoma/store/old_Varscan2/'
script='/work/szlab/kh31516_Lab_Share_script'
Total_file_name='/scratch/kh31516/Original_Melanoma/source/combined_melanoma.txt'

#result='/scratch/kh31516/Pan_cancer/glioma/results/WGS_WES/finished/${line1}'
#source_data='/scratch/kh31516/Original_Mammary/results/${line1}'

### filter Lofreq #####
ml Anaconda3/2019.07
cd ${Lofreq_folder}
### Limit the indels in the CDS region ###

while read line1 line2 line3;
do
    python $script/lofreq_limitCDS.py \
    $Lofreq_folder/${line1}/${line1}_lofreq.vcfsomatic_final_minus-dbsnp.snvs.vcf.gz \
    $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
    $Lofreq_folder/${line1}/${line1}_Lofreq_SNV.vcf_canFam3.1.99_CDS.vcf

    python $script/lofreq_limitCDS.py \
    $Lofreq_folder/${line1}/${line1}_lofreq.vcfsomatic_final_minus-dbsnp.indels.vcf.gz \
    $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
    $Lofreq_folder/${line1}/${line1}_Lofreq_indel.vcf_canFam3.1.99_CDS.vcf
done < $Total_file_name




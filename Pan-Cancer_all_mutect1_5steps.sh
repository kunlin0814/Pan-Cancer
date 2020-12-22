#! /bin/bash

: "
while read line1 line2 line3
do 
    mkdir /scratch/kh31516/Original_Melanoma/store/PRJEB12081/VAF/new_VAF/$line1
    cd /scratch/kh31516/Original_Melanoma/store/PRJEB12081/VAF/Mutect/$line1
    cat *_after.txt | cut -f1-3 > /scratch/kh31516/Original_Melanoma/store/PRJEB12081/VAF/new_VAF/$line1/${line1}_vaf_after.txt
    cat *_before.txt | cut -f1-3 > /scratch/kh31516/Original_Melanoma/store/PRJEB12081/VAF/new_VAF/$line1/${line1}_vaf_before.txt

done < /scratch/kh31516/Original_Melanoma/source/combined_melanoma.txt 
"

ml Anaconda3/2020.02
## 5 steps filtering

scripts='/home/kh31516/kh31516_Lab_Share_script'

#MC
total_file='/scratch/kh31516/Original_Mammary/source/K9_WXS_mammary_SRR_pairs.txt'
mutect_folder='/scratch/kh31516/Original_Mammary/store/PRJNA489159/newMutectSNP'
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/MT_VAF/Mutect1'

while read line1 line2 line3
do
    #echo $line1
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file

#OM
total_file='/scratch/kh31516/Original_Melanoma/source/combined_melanoma.txt'
mutect_folder='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/newMutectSNP'
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/OM_VAF/Mutect1'


while read line1 line2 line3
do
   mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


# GLM
total_file='/scratch/kh31516/Pan_cancer/glioma/source/Glioma_WES_Pair.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/glioma/results/store/WES/new_Mutect_SNP'
vaf_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WES/Mutect1VAF'


while read line1 line2 line3
do
     mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file



# HSA
total_file='/scratch/kh31516/Pan_cancer/HSA/source/PRJNA552034-HSA-Normal-Tumor-pair'
mutect_folder='/scratch/tw71066/All_NEW_WES/PRJNA552034-HSA/Mutation/new_mutect_results'
vaf_out='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA552034/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
done < $total_file



# LYM
total_file='/scratch/kh31516/Pan_cancer/lymphoma/source/Lym_pairs.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/lymphoma/store/Mutect_Mar25/'
vaf_out='/scratch/kh31516/Pan_cancer/lymphoma/store/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


#OSA
total_file='/scratch/kh31516/Pan_cancer/Osteo/source/Bur_osteo_pair'
mutect_folder='/scratch/kh31516/Pan_cancer/Osteo/store/PRJNA391455/Mutect_Mar25_results'
vaf_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJNA391455/Mutect1VAF'



while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


#UCL

total_file='/scratch/kh31516/Pan_cancer/others/source/Unclassified_pairs'
mutect_folder='/scratch/kh31516/Pan_cancer/others/store/Mutect_Mar25'
vaf_out='/scratch/kh31516/Pan_cancer/others/store/Mutect1VAF'



while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file
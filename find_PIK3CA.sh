#!/bin/bash
#PBS -q batch
#PBS -N new-bwa-java-${line1}
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=55gb

results='/scratch/kh31516/Original_Mammary/store/PRJNA552905/Mutect/'

while read line1 line2 line3 ;
do
    check=$(cat $results/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function_WithGeneName | grep "PIK3CA")
    if [ ! -z "$check" ]; then
        
        gene=$(cat $results/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function_WithGeneName | grep "PIK3CA" | awk -F "\t" '{print $4}' | awk -F ":" '{print $5}')

        printf "%s\t%s\n" ${line1} ${gene} >> $results/PIK3CA_location.txt
        #printf "%s\t%s\n" "${line1} "${gene}" 
    fi    
       

done < /scratch/kh31516/Original_Mammary/source/PRJNA552905_MC_pairs.txt

 #>> $results/PIK3CA_location.txt

while read line1 line2 line3 ;
do
    results='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA417727/Mutect'
    check=$(cat $results/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function_WithGeneName | grep "TP53")
    if [ ! -z "$check" ]; then
        gene=$(cat $results/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function_WithGeneName | grep "TP53" | awk -F "\t" '{print $4}' | awk -F ":" '{print $5}')

        printf "%s\t%s\n" ${line1} ${gene} >> $results/HSA_TP53_location.txt

    fi
done <  /scratch/kh31516/Pan_cancer/HSA/source/PRJNA417727_HSA_WES.txt

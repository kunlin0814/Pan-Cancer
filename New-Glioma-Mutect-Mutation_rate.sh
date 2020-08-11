#!/bin/bash
#PBS -N retro_non_retro-glioma
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=30gb
#PBS -q batch


source='/work/szlab/Retro_gene'
result='/scratch/kh31516/Pan_cancer/glioma_new_mut/glioma'
Total_file_name='/scratch/kh31516/Pan_cancer/source/glioma_cases'

module load Python/3.7.0-foss-2018a
#module load Anaconda3/2018.12
#source activate py35

cd  $result/

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "file_name" "retro_mutation_rate" "retor_PASS" "retro_callable" "non_retro_mutation_rate" \
"non_retro_PASS" "nonr_retro_callable" "combine_mutation_rate" "combine_PASS" "combine_callable" >> $result/New_Mut_retro_non_retro_melanoma_summary.txt
while read line;
do
    cd  $result/${line}
    #awk '$7 == "PASS" {print $0}' ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut
    printf "%s\t%s\n" "${line}" "$result/${line}/${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt" > $result/${line}/sample_wigs_file.txt
    a=$(cat ${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1" | wc -l | bc) # denominator, all callable
    b=$(cat ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut | grep -v "^#" | grep -i "PASS"| wc -l |bc) #numerator, all PASS
    c=$(python $source/find_vcf.py $source/new_retro_gene_list_information_CanFam3.1.99gtf.txt $result/${line}/${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut) # retro numerator
    java -cp $source/ GetCallableCounts $result/${line}/sample_wigs_file.txt $result/${line}/retro_callable_summary.txt # give the retro gene callable summary
    d=$(python $source/calcaulte_retro_callable.py $result/${line}/retro_callable_summary.txt) # retro denominator
    e=$((b-c)) # non-retro numerator
    f=$((a-d)) # non-retro denominator
    g=$(echo "$((c))*1000000/$((d))" | bc -l) # retro mutation rate per million
    h=$(echo "$((e))*1000000/$((f))" | bc -l) # non-retro mutation rate per million
    i=$(echo "$((b))*1000000/$((a))" | bc -l)
    printf  "%s\t%4f\t%d\t%d\t%4f\t%d\t%d\t%4f\t%d\t%d\n" "${line}" "${g}" "${c}" "${d}" "${h}" "${e}" "${f}" \
    "${i}" "${b}" "${a}" >> $result/New_Mut_retro_non_retro_melanoma_summary.txt
    rm $result/${line}/sample_wigs_file.txt
    rm $result/${line}/retro_callable_summary.txt
done < $Total_file_name
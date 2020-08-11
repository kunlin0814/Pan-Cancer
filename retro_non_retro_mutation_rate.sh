#!/bin/bash
#PBS -N retro_non_retro
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=30gb
#PBS -q batch


module load Python/3.7.0-foss-2018a

cd  /scratch/kh31516/Original_Melanoma/

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "file_name" "retro_mutation_rate" "retor_PASS" "retro_callable" "non_retro_mutation_rate" \
"non_retro_PASS" "nonr_retro_callable" "combine_mutation_rate" "combine_PASS" "combine_callable" >> /scratch/kh31516/Original_Melanoma/retro_non_retro_melanoma_summary.txt
while read line;
do
    cd  /scratch/kh31516/Original_Melanoma/result/${line}
    awk '$7 == "PASS" {print $0}' ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
    printf "%s\t%s\n" "${line}" "/scratch/kh31516/Original_Melanoma/result/${line}/${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt" > /scratch/kh31516/Original_Melanoma/result/${line}/sample_wigs_file.txt
    a=$(cat ${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1" | wc -l | bc) # denominator, all callable
    b=$(cat ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf | grep -v "^#" | grep -i "PASS"| wc -l |bc) #numerator, all PASS
    c=$(python3.7 /scratch/kh31516/Original_Mammary/analysis_scripts/find_vcf.py /scratch/kh31516/Original_Mammary/source/retro_gene_list_information.txt /scratch/kh31516/Original_Melanoma/result/${line}/${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS) # retro numerator
    java -cp /scratch/kh31516/Original_Mammary/analysis_scripts/ GetCallableCounts /scratch/kh31516/Original_Melanoma/result/${line}/sample_wigs_file.txt /scratch/kh31516/Original_Melanoma/result/${line}/retro_callable_summary.txt # give the retro gene callable summary
    d=$(python3.7 /scratch/kh31516/Original_Mammary/analysis_scripts/calcaulte_retro_callable.py /scratch/kh31516/Original_Melanoma/result/${line}/retro_callable_summary.txt) # retro denominator
    e=$((b-c)) # non-retro numerator
    f=$((a-d)) # non-retro denominator
    g=$(echo "$((c))*1000000/$((d))" | bc -l) # retro mutation rate per million
    h=$(echo "$((e))*1000000/$((f))" | bc -l) # non-retro mutation rate per million
    i=$(echo "$((b))*1000000/$((a))" | bc -l)
    printf  "%s\t%4f\t%d\t%d\t%4f\t%d\t%d\t%4f\t%d\t%d\n" "${line}" "${g}" "${c}" "${d}" "${h}" "${e}" "${f}" \
    "${i}" "${b}" "${a}" >> /scratch/kh31516/Original_Melanoma/retro_non_retro_melanoma_summary.txt

done < /scratch/kh31516/Original_Melanoma/source/melanoma_cases.txt
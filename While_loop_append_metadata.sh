#! /bin/bash

while read line ; 
do
cd /scratch/kh31516/Mammary/results/${line}
python /scratch/kh31516/Mammary/scripts/Append_annovar_file_with_MetaInfo.py ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName /scratch/kh31516/Mammary/source/Mammary_WES.txt $line
done < /scratch/kh31516/Mammary/source/Mammary_cases.txt

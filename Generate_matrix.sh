
script='/scratch/kh31516/Mammary/scripts/DNA_repair'
data='/scratch/kh31516/Mammary/results/'

cd $script
python $script/Annovar_DNArepair_summary.py $data/CMT-100/CMT-100_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName name CMT-100 > /scratch/kh31516/Mammary/results/Mammary_mutect_DNArepair_somatic.txt

while read line
do
python $script/Annovar_DNArepair_summary.py $data/$line/$line'_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName' mut $line >> /scratch/kh31516/Mammary/results/Mammary_mutect_DNArepair_somatic.txt
done < /scratch/kh31516/Mammary/source/Mammary_cases.txt



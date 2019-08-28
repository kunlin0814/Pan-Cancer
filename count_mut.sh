#PBS -S /bin/bash
#PBS -q batch
#PBS -N count_mut_gene
#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00
#PBS -l mem=40gb
#PBS -M kh31516@uga.edu 
#PBS -m ae

#for i in /scratch/kh31516/Mammary/results/CMT*/CMT-*_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName
while read line;
do 
    cd /scratch/kh31516/Mammary/results/$line
    python /scratch/kh31516/Mammary/CRC_find_mut.py non_syn_${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName /scratch/kh31516/Mammary/total_CRC_pathway.txt
done < /scratch/kh31516/Mammary/source/Mammary_cases.txt

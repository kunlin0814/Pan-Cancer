#!/bin/bash
#SBATCH --job-name=lab_5_Summarized_BWA_WES         # Job name (Summarized_BWA_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=55G                   # Job memory limit (10 GB)
#SBATCH --time=10:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=lab_5_Summarized_BWA_WES.%j.out    # Standard output log
#SBATCH --error=lab_5_Summarized_BWA_WES.%j.err     # Standard error log


# OM

total_file='/scratch/kh31516/Original_Melanoma/source/combined_melanoma.txt'
output_file='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/totalSangerOM_mutation.txt'
while read line1 line2 line3
do
    mutect_folder='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/newMutectSNP'/$line1
    
    cd $mutect_folder

    cat ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut | awk -F "\t" '{print $1,$2,$4,$5}' | \
    sed 's/ /"\t"/g' |sed 's/"//g' | awk -v awkvar="$line1" -F"\t" 'BEGIN { OFS = "\t" } {$5=awkvar; print}' >> $output_file
    
    strelka_folder='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/strelka_result/'${line1}-demo_somatic/results/variants

    indelFile=${strelka_folder}/somatic.indels.vcf_canFam3.1.99_CDS-PASS
    
    if [ -s $indelFile ]; then
        cat $indelFile |  awk -F "\t" '{print $1,$2,$4,$5}' | \
        sed 's/ /"\t"/g' |sed 's/"//g' | awk -v awkvar="$line1" -F"\t" 'BEGIN { OFS = "\t" } {$5=awkvar; print}' >> $output_file
    fi

done < $total_file


# MT
total_file='/scratch/kh31516/Original_Mammary/source/K9_WXS_mammary_SRR_pairs.txt'
output_file='/scratch/kh31516/Original_Mammary/store/PRJNA489159/totalMT_mutation.txt'
while read line1 line2 line3
do
    mutect_folder='/scratch/kh31516/Original_Mammary/store/PRJNA489159/newMutectSNP'/$line1
    
    cd $mutect_folder

    cat ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut | awk -F "\t" '{print $1,$2,$4,$5}' | \
    sed 's/ /"\t"/g' |sed 's/"//g' | awk -v awkvar="$line1" -F"\t" 'BEGIN { OFS = "\t" } {$5=awkvar; print}' >> $output_file
    
    strelka_folder='/scratch/kh31516/Original_Mammary/store/PRJNA489159/strelka_result'/${line1}-demo_somatic/results/variants

    indelFile=${strelka_folder}/somatic.indels.vcf_canFam3.1.99_CDS-PASS
    
    if [ -s $indelFile ]; then
        cat $indelFile |  awk -F "\t" '{print $1,$2,$4,$5}' | \
        sed 's/ /"\t"/g' |sed 's/"//g' | awk -v awkvar="$line1" -F"\t" 'BEGIN { OFS = "\t" } {$5=awkvar; print}' >> $output_file
    fi

done < $total_file
#!/bin/bash
#SBATCH --job-name=MT_ensembl_callable         # Job name (MT_ensembl_callable)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=50G                   # Job memory limit (10 GB)
#SBATCH --time=48:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=MT_ensembl_callable.%j.out    # Standard output log
#SBATCH --error=MT_ensembl_callable.%j.err     # Standard error log


tumor_type="GLM"
scripts="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_script"
total_file='/scratch/kh31516/Pan_cancer/glioma/source/Glioma_WES_Pair.txt'
# mutect results from Burair_filtering3
mutect_folder="/scratch/kh31516/Pan_cancer/Burair_filtering3/GLM/Mutect"
ensembl_source="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_source/Pan_cancer_sign_target_ensembl_location.txt"
wig_callable="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_source/OM_wig_samples.txt"
Retro_gene_source='/home/kh31516/kh31516_Lab_Share_script/Retro_gene_source'
Retro_gene_script='/home/kh31516/kh31516_Lab_Share_script/Retro_gene_script'
output_folder='/scratch/kh31516/Pan_cancer/Burair_filtering3/GLM/Mutect'
wig_folder="/scratch/kh31516/Pan_cancer/glioma/results/store/WES/new_Mutect_SNP"
total_mutation_number_file="/scratch/kh31516/total_${tumor_type}_mutation_inlist.txt"
total_genome_mutation_number_file="/scratch/kh31516/total_${tumor_type}_mutation_genomewide.txt"
ml Anaconda3/2020.02

while read line line1 line2;
do
    echo "processing $line"

    cd  $mutect_folder/$line

    # candidate_info = sys.argv[1]
    # mut_gene_info = sys.argv[2]
    # sample_name = sys.argv[3]
    # output = open(sys.argv[4], 'w')

    python $scripts/find_ensmbl_mutation.py \
    $ensembl_source \
    $line.Mutect.vcf-5steps_orientBias3-avinput.exonic_variant_function_WithGeneName \
    $line \
    $line_filtering_mutation_number_inlist.txt

    
    cat $line_filtering_mutation_number_inlist.txt >> $total_mutation_number_file


done < $total_file


ml Java

        # String annotationFile = args[0];
        # //"C:/Users/abc73/Desktop/TestTMB/Canis_familiaris.CanFam3.1.99.chr.gtf";
        # String geneListFile = args[1];
        # //"C:/Users/abc73/Desktop/TestTMB/Pan_cancer_sign_target_ensembl_location.txt";
        # String sampleWigFile = args[2];
        # //"C:/Users/abc73/Desktop/TestTMB/sample_wig.txt";
        # String outputFile = args[3];

total_mutation_callable_file="/scratch/kh31516/total_${tumor_type}_callable_mutation_inlist.txt"

java -Xmx32g -cp $scripts/ GetCallableCounts /work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.chr.gtf \
$ensembl_source \
$wig_callable $total_mutation_callable_file

gzip $total_mutation_number_file
gzip $total_mutation_callable_file

cd  ${mutect_folder}

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "file_name" "retro_mutation_rate" "retor_PASS" "retro_callable" "non_retro_mutation_rate" \
"non_retro_PASS" "nonr_retro_callable" "combine_mutation_rate" "combine_PASS" "combine_callable" "cancer_type" >> ${total_genome_mutation_number_file}

while read line line1 line2;
do
    mkdir -p ${output_folder}/${line}
    cd ${wig_folder}/${line}
    #awk '$7 == "PASS" {print $0}' ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${line}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut
    printf "%s\t%s\n" "${line}" "${wig_folder}/${line}/${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt" > ${output_folder}/${line}/sample_wigs_file.txt
    a=$(cat ${line}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1" | wc -l | bc) # denominator, all callable
    cd ${mutect_folder}/${line}/
    b=$(cat ${line}.Mutect.vcf-5steps_orientBias3 | grep -v "^#" | grep -i "PASS"| wc -l |bc) #numerator, all PASS
    c=$(python $Retro_gene_script/find_vcf.py $Retro_gene_source/new_retro_gene_list_information_CanFam3.1.99gtf.txt ${mutect_folder}/${line}/${line}.Mutect.vcf-5steps_orientBias3) # retro numerator
    
    java -cp $scripts/ GetCallableCounts /work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.chr.gtf \
    $Retro_gene_source/new_retro_gene_list_information_CanFam3.1.99gtf.txt \
    ${output_folder}/${line}/sample_wigs_file.txt \
    ${output_folder}/${line}/retro_callable_summary.txt # give the retro gene callable summary
    
    d=$(python $Retro_gene_script/calcaulte_retro_callable.py ${output_folder}/${line}/retro_callable_summary.txt) # retro denominator
    e=$((b-c)) # non-retro numerator
    f=$((a-d)) # non-retro denominator
    g=$(echo "$((c))*1000000/$((d))" | bc -l) # retro mutation rate per million
    h=$(echo "$((e))*1000000/$((f))" | bc -l) # non-retro mutation rate per million
    i=$(echo "$((b))*1000000/$((a))" | bc -l)
    printf  "%s\t%4f\t%d\t%d\t%4f\t%d\t%d\t%4f\t%d\t%d\t%s\n" "${line}" "${g}" "${c}" "${d}" "${h}" "${e}" "${f}" \
    "${i}" "${b}" "${a}" "${tumor_type}" >> $total_genome_mutation_number_file
    
    rm ${output_folder}/${line}/sample_wigs_file.txt
    #rm -r ${output_folder}/${line}/retro_callable_summary.txt
done < $total_file

gzip $total_genome_mutation_number_file
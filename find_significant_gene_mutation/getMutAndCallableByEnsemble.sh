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


tumor_type="MT"
scripts="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_script"
total_file='/scratch/kh31516/Original_Mammary/source/K9_WXS_mammary_SRR_pairs.txt'
# mutect results from Burair_filtering3
mutect_folder="/scratch/kh31516/Pan_cancer/Burair_filtering3/MT_Korean/Mutect"
ensembl_source="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_source/Pan_cancer_sign_target_ensembl_location.txt"
wig_callable="/home/kh31516/kh31516_Lab_Share_script/Retro_gene_source/MT_wig_samples.txt"
total_mutation_number_file="/scratch/kh31516/total_${tumor_type}_mutation_inlist.txt"

ml Anaconda3/2020.02

# MT Ko
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


## MT SUN

total_file='/scratch/kh31516/Original_Mammary/source/PRJNA552905_MC_pairs.txt'
# mutect results from Burair_filtering3
mutect_folder="/scratch/kh31516/Pan_cancer/Burair_filtering3/MT_SNU/Mutect"


while read line line1 line2;
do
    echo "processing $line"
    cd  $mutect_folder/$line

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

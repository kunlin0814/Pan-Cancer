#!/bin/bash
#SBATCH --job-name=filtering_OM_PUB         # Job name (filtering_OM_PUB)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=40G                   # Job memory limit (10 GB)
#SBATCH --time=90:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=filtering_OM_PUB.%j.out    # Standard output log
#SBATCH --error=filtering_OM_PUB.%j.err     # Standard error log

## The script take CDS mutation file from OM, 
## The OM mutation first filtering with only CDS region and then filtering with DBSNP

reference='/work/szlab/Lab_shared_PanCancer/source'
file_folder='/scratch/kh31516' 
script='/home/kh31516/kh31516_Lab_Share_script/Mutect2'
output_folder='/scratch/kh31516'
#result='/scratch/kh31516/Pan_cancer/glioma/results/WGS_WES/finished/CMT-2'
#source_data='/scratch/kh31516/Original_Mammary/results/CMT-2'

ml Java

cd ${file_folder}
### filter  #####
## with DbSNP

java -Xmx32g -cp  $script/ DbSNP_publish_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${file_folder}/sanger_mut_file_before_DBSNP_03_23.txt \
${file_folder}/DbSNP_sanger_mut_file_after_DBSNP_03_23.txt

## with PON

java -Xmx32g -cp  $script/ DbSNP_publish_filtering \
$reference/pon.vcf.txt \
${file_folder}/DbSNP_sanger_mut_file_after_DBSNP_03_23.txt \
${file_folder}/PON_DbSNP_sanger_mut_file_after_DBSNP_03_23.txt


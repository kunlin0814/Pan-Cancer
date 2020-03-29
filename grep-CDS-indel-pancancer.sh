#!/bin/bash
#PBS -N grep-CDS-indel-pancancer
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -q batch

 
cd /scratch/kh31516/Pan_cancer/
: "
cp -r /scratch/baa20876/workDir/pancancer/WES/PRJNA247493/lymphoma/results/strelka/* /scratch/kh31516/Pan_cancer/other_Strelka/Lym/
cp -r /scratch/baa20876/workDir/pancancer/WES/PRJNA247493/other/results/strelka/* /scratch/kh31516/Pan_cancer/other_Strelka/Others/
cp -r /scratch/baa20876/workDir/pancancer/WES/PRJNA247493/osteosarcoma/strelka/* /scratch/kh31516/Pan_cancer/other_Strelka/Osteo
"

### Bur_Osteo ###
source_file='/scratch/kh31516/Pan_cancer/Bur_Osteo/source/Bur_Osteo_cases.txt'

cd /scratch/kh31516/Pan_cancer/other_Strelka/Osteo

python /scratch/kh31516/Limit_vcf_to_CDS.py \
"/scratch/kh31516/Pan_cancer/other_Strelka/Osteo/*/results/variants/somatic.indels.vcf.gz" \
/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list

while read line;
do
    cd /scratch/kh31516/Pan_cancer/other_Strelka/Osteo/${line}/results/variants
    num_indel=$(cat somatic.indels.vcf_canFam3.1.99_CDS | grep -v "^##" | grep -w "PASS" | wc -l)
    printf "%s\t%d\n" "${line}" "${num_indel}" >> /scratch/kh31516/Pan_cancer/other_Strelka/Bur_Osteo-CDS-indel_mutation_number.txt
done < ${source_file}

### Bur_Lym ###
cd /scratch/kh31516/Pan_cancer/other_Strelka/Lym/
python /scratch/kh31516/Limit_vcf_to_CDS.py \
"/scratch/kh31516/Pan_cancer/other_Strelka/Lym/*/results/variants/somatic.indels.vcf.gz" \
/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list

source_file='/scratch/kh31516/Pan_cancer/source/lym_cases'

cd /scratch/kh31516/Pan_cancer/other_Strelka/Lym/

while read line;
do
    cd /scratch/kh31516/Pan_cancer/other_Strelka/Lym/${line}/results/variants
    num_indel=$(cat somatic.indels.vcf_canFam3.1.99_CDS | grep -v "^##" | grep -w "PASS" | wc -l)
    printf "%s\t%d\n" "${line}" "${num_indel}" >> /scratch/kh31516/Pan_cancer/other_Strelka/Lymphoma-CDS-indel_mutation_number.txt
done < ${source_file}


### Bur_Others ###
cd /scratch/kh31516/Pan_cancer/other_Strelka/Others/
python /scratch/kh31516/Limit_vcf_to_CDS.py \
"/scratch/kh31516/Pan_cancer/other_Strelka/Others/*/results/variants/somatic.indels.vcf.gz" \
/work/szlab/Lab_shared_PanCancer/source/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list

source_file='/scratch/kh31516/Pan_cancer/source/other_cases'

cd /scratch/kh31516/Pan_cancer/other_Strelka/Others/

while read line;
do
    cd /scratch/kh31516/Pan_cancer/other_Strelka/Others/${line}/results/variants
    num_indel=$(cat somatic.indels.vcf_canFam3.1.99_CDS | grep -v "^##" | grep -w "PASS" | wc -l)
    printf "%s\t%d\n" "${line}" "${num_indel}" >> /scratch/kh31516/Pan_cancer/other_Strelka/Others-CDS-indel_mutation_number.txt
done < ${source_file}
#!/bin/bash
#SBATCH --job-name=i_03A6_Germline_Callable         # Job name (i_03A6)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU co4es per task (4)
#SBATCH --mem=70G                   # Job memory limit (10 GB)
#SBATCH --time=160:00:00              # Time 1it hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=i_03A6.%j.out    # Standard output log
#SBATCH --error=i_03A6.%j.err     # Standard1ror4log



sample_name='i_03A6'
Normal_Run='SRR10362478'
Tumor_Run='SRR10362479'
results='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/results/'${sample_name} #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/'${sample_name}/${Normal_Run}
dataT='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/data/'${sample_name}/${Tumor_Run}
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/work/szlab/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect/'${sample_name}
MuTect2_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Mutect2/'${sample_name}
Germline_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/Germline/'${sample_name}
DepthOfCoverage='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/DepthOfCoverage/'${sample_name}
strelka_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WGS/strelka/'


############################## Germline mutation preparation Start ######################################
# GATK 3 
ml GATK/3.8-1-Java-1.8.0_144

cd $Germline_out

# Variant calling
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/canFam3.fa \
-V ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

################ Annovar ###################

cd ${Germline_out}
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

awk '$7 == "PASS" {print $0}' ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS

# annovar input preparation
module load Perl/5.26.1-GCCcore-6.4.0

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

perl $annovar_index/convert2annovar.pl -format vcf4old ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS \
> ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

perl $annovar_index/annotate_variation.pl --buildver canFam3 ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput $annovar_index

# add gene names to annovar output
ml Anaconda3/2020.02
cd ${Germline_out}
python2 $script/Add_GeneName_N_Signature.py ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

python2 $script/Add_GeneName_N_Signature.py ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput.exonic_variant_function \
${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


#### GATK callable ####
cd ${Germline_out}
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Normal_Run}_table.txt \
-o ${Germline_out}/${Normal_Run}_callable_status.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Tumor_Run}_table.txt \
-o ${Germline_out}/${Tumor_Run}_callable_status.bed 

############################## Germline mutation preparation End ######################################


#### GATK DepthofCoverage ####
cd ${DepthOfCoverage}
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.interval_list \
-o ${DepthOfCoverage}/${Normal_Run}_DepthofCoverage_CDS.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.interval_list \
-o ${DepthOfCoverage}/${Tumor_Run}_DepthofCoverage_CDS.bed 


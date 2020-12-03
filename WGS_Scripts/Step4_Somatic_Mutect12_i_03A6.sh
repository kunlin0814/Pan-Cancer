#!/bin/bash
#SBATCH --job-name=i_03A6_Mutect12         # Job name (i_03A6)
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


######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 

cd ${MuTect_out}
module load MuTect/1.1.7-Java-1.7.0_80

java -jar $EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.interval_list \
--input_file:normal ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf


awk '$7 == "PASS" {print $0}' ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS


ml Anaconda3/2020.02

python $script/Filter_MutectStat_5steps.py \
${MuTect_out}/${sample_name}_PASS.stat \
${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
${MuTect_out}/${sample_name}_vaf_before.txt \
${MuTect_out}/${sample_name}_vaf_after.txt \
${MuTect_out}/${sample_name}_whyout.txt 


## No 5 steps filtering  

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# add gene name
python2 $script/Add_GeneName_N_Signature.py ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


## With 5 steps filtering 

### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${MuTect_out}/${sample_name}_PASS.stat

## 5 steps filtering created XXX_filterMut file 

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python2 $script/Add_GeneName_N_Signature.py ${MuTect_out}/${sample_name}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


####### Mutect2 #######

cd ${MuTect2_out}

ml SAMtools/1.9-GCC-8.3.0

### Change Bam file Read groups ###
## if you need to re-header and re-index , the output re-alignment file must be in the different location and can't overwrite
##
samtools view -H ${results}/SRR10362478_rg_added_sorted_dedupped_removed.realigned.bam | sed -e "s/SM:20/SM:SRR10362478_normal/" | \
samtools reheader - ${results}/SRR10362478_rg_added_sorted_dedupped_removed.realigned.bam > ${MuTect2_out}/SRR10362478_rg_added_sorted_dedupped_removed.realigned.bam

samtools view -H ${results}/SRR10362479_rg_added_sorted_dedupped_removed.realigned.bam | sed -e "s/SM:20/SM:SRR10362479_tumor/" | \
samtools reheader - ${results}/SRR10362479_rg_added_sorted_dedupped_removed.realigned.bam > ${MuTect2_out}/SRR10362479_rg_added_sorted_dedupped_removed.realigned.bam


mv ${MuTect2_out}/SRR10362478_rg_added_sorted_dedupped_removed.realigned.bam ${results}
mv ${MuTect2_out}/SRR10362479_rg_added_sorted_dedupped_removed.realigned.bam ${results}

## GATK re index 
module load picard/2.16.0-Java-1.8.0_144

java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
I=${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam

java -jar $EBROOTPICARD/picard.jar BuildBamIndex \
I=${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

cd ${MuTect2_out}
ml SAMtools/1.9-GCC-8.3.0
ml GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

Normal_sample=$(samtools view -H ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--callable-depth 8 \
--dont-use-soft-clipped-bases true \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${MuTect2_out}/${sample_name}-f1r2.tar.gz \
-L $MuTect2_source/Canis_familiaris.CanFam3.interval.list \
-O ${MuTect2_out}/${sample_name}_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/${sample_name}-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/${sample_name}_MuTect2_GATK4_noDBSNP.vcf \
-O ${MuTect2_out}/filtered-${sample_name}_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-${sample_name}_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf


## with PON
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf

### Mutect2 5 steps filtering

cd ${MuTect2_out}

ml Anaconda3/2020.02

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS

awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS \
${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS \
${MuTect2_out}/${sample_name}_VAF_Before.txt \
${MuTect2_out}/${sample_name}_VAF_After.txt \
${MuTect2_out}/${sample_name}_whyout.txt


####### Annovar #######
module load Perl/5.26.1-GCCcore-6.4.0

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS_filteredMut > ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python2 $script/Add_GeneName_N_Signature.py ${MuTect2_out}/PON_DbSNP_filtered-${sample_name}_MuTect2_GATK4.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



## Strelka
cd ${results}
module load SAMtools/1.9-GCC-8.3.0 
#samtools index ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam
#samtools index ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

${script}/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--tumorBam ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--referenceFasta ${reference}/canFam3.fa \
--runDir $strelka_out/${sample_name}

# your result will be generated in demo_somatic/results/variants
$strelka_out/${sample_name}/runWorkflow.py -m local -j 20

: "
## limit strelka result into CDS region
## WGS current don't consdier CDS

python2 $script/Limit_vcf_to_CDS.py $strelka_out/results/variants/somatic.indels.vcf.gz ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list

"

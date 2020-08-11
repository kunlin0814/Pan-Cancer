#!/bin/bash
#PBS -q batch
#PBS -N i_3F8C-Mutect2-Varscan2-somatic_seq-filtering
#PBS -l nodes=1:ppn=4
#PBS -l walltime=80:00:00
#PBS -l mem=60gb


MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
reference='/work/szlab/Lab_shared_PanCancer/source'
Mutect2_result='/scratch/kh31516/Pan_cancer/glioma/results/WES_Mutect2/i_3F8C' 
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
Lofreq_folder='/scratch/kh31516/Pan_cancer/glioma/results/store/lofreq/i_3F8C'
VarScan2_folder='/scratch/kh31516/Pan_cancer/glioma/results/store/Varscan2/i_3F8C'
somatic_seq_output_folder='/scratch/kh31516/Pan_cancer/glioma/results/store/somatic_seq/i_3F8C'
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
script='/work/szlab/kh31516_Lab_Share_script'
Glioma_WES_result='/scratch/kh31516/Pan_cancer/glioma/results/WES/i_3F8C'
#source_data='/scratch/kh31516/Original_Mammary/results/i_3F8C'


mkdir -p ${Mutect2_result}
cd ${Mutect2_result}

ml SAMtools/1.9-GCC-8.2.0-2.31.1

### Change Bam file Read groups ###
samtools view -H ${Glioma_WES_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam | sed -e "s/SM:20/SM:SRR10351806_normal/" | \
samtools reheader - ${Glioma_WES_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam > ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam

samtools view -H ${Glioma_WES_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam | sed -e "s/SM:20/SM:SRR10351807_tumor/" | \
samtools reheader - ${Glioma_WES_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam > ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam


### Create an Index for GATK use ###
cd ${Mutect2_result}
module load picard/2.16.0-Java-1.8.0_144

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex \
I=${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex \
I=${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam


Normal_sample=$(samtools view -H ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)


### LoFreq ###
mkdir -p $Lofreq_folder
cd $Lofreq_folder
ml bzip2/1.0.6-foss-2016b
module load LoFreq/2.1.2-foss-2016b-Python-2.7.14

lofreq somatic -n ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam \
-t ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam \
-f $reference/canFam3.fa \
--threads 4 -o $Lofreq_folder/i_03A6_lofreq.vcf \
-d $Lofeq_source/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf.gz \
--call-indels \
--min-cov 7 \
--verbose



## change the GATK version and Run Mutect2 ##
cd ${Mutect2_result}
module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam \
--callable-depth 8 \
-normal $Normal_sample \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${Mutect2_result}/i_3F8C-f1r2.tar.gz \
-L $GATK4/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${Mutect2_result}/i_3F8C_MuTect2_GATK4_noDBSNP.vcf

## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${Mutect2_result}/i_3F8C-f1r2.tar.gz -O ${Mutect2_result}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${Mutect2_result}/read-orientation-model.tar.gz \
-V ${Mutect2_result}/i_3F8C_MuTect2_GATK4_noDBSNP.vcf \
-O ${Mutect2_result}/filtered-i_3F8C_MuTect2_GATK4_noDBSNP.vcf


######### Varscan2 ###########
ml SAMtools/1.9-GCC-8.2.0-2.31.1
ml VarScan/2.4.2-Java-1.8.0_144

mkdir -p $VarScan2_folder
cd ${VarScan2_folder}

samtools mpileup -f $reference/canFam3.fa -q 1 -B ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam \
${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam > ${VarScan2_folder}/SRR10351806-SRR10351807.mpileup

time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somatic \
${VarScan2_folder}/SRR10351806-SRR10351807.mpileup \
${VarScan2_folder}/1mpileup-i_3F8C_removed_realigned_Varscan \
--mpileup 1 \
--min-coverage 8 \
--min-coverage-normal 6 \
--min-coverage-tumor 8 \
-–min-reads2 2 \
-–min-avg-qual 15 \
-–min-var-freq 0.08 \
-–min-freq-for-hom 0.75 \
–-tumor-purity 1.0 \
-–strand-filter 1 \
--–somatic-p-value 0.05 \
--output-vcf 1

#### Filter the somatic mutation derived from Varscan2 ####
time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter ${VarScan2_folder}/1mpileup-i_3F8C_removed_realigned_Varscan.snp.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file ${VarScan2_folder}/i_3F8C_Varscan2.vcf.snp_filtered.vcf

time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter ${VarScan2_folder}/1mpileup-i_3F8C_removed_realigned_Varscan.indel.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file ${VarScan2_folder}/i_3F8C_Varscan2.vcf.indel_filtered.vcf


############################### Somatic Seq #################################
### Somatic Seq , needs the filtered Mutect2, filtered lofreq and filtered varscan2 vcf file, 
## The order of 3 VCF files cannot change, The orders are Mutect2 , lofreq, varscan2
## ( it will affect the further the analysis prioritized by Mutect2) ###

mkdir -p $somatic_seq_output_folder
cd  $somatic_seq_output_folder
module load somaticseq/3.4.1_conda
module load BEDTools/2.28.0-foss-2018a
source activate ${SOMATICSEQROOT}

## all the input vcf were filtered before somaticseq, 
# these input files didn't not go through dbsnp nor limit to CDS region

somaticseq_parallel.py \
--output-directory ${somatic_seq_output_folder} \
--genome-reference ${reference}/canFam3.fa \
--inclusion-region ${reference}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_tab.bed \
--algorithm ada \
--threads 4 \
paired \
--tumor-bam-file    ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam \
--normal-bam-file   ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam \
--mutect2-vcf       ${Mutect2_result}/filtered-i_3F8C_MuTect2_GATK4_noDBSNP.vcf \
--lofreq-snv        ${Lofreq_folder}/i_3F8C_lofreq.vcfsomatic_final.snvs.vcf.gz \
--lofreq-indel      ${Lofreq_folder}/i_3F8C_lofreq.vcfsomatic_final.indels.vcf.gz \
--varscan-snv       ${VarScan2_folder}/i_3F8C_Varscan2.vcf.snp_filtered.vcf \
--varscan-indel     ${VarScan2_folder}/i_3F8C_Varscan2.vcf.indel_filtered.vcf 


####### filtering the consensus VCF files created by somaticseq with Dbsnp #########

cd $somatic_seq_output_folder

java -cp  $script/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${somatic_seq_output_folder}/Consensus.sSNV.vcf \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sSNV.vcf

java -cp  $script/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${somatic_seq_output_folder}/Consensus.sINDEL.vcf \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sINDEL.vcf

### Change the filter status and prioritize by Mutect2 

cd $somatic_seq_output_folder

python ${script}/prioritizedByMutect2.py \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sSNV.vcf \
${somatic_seq_output_folder}/Mutect2_Prioritized_DbSNPFiltering_Consensus.sSNV.vcf

python ${script}/prioritizedByMutect2.py \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sINDEL.vcf \
${somatic_seq_output_folder}/Mutect2_Prioritized_DbSNPFiltering_Consensus.sINDEL.vcf

## remove the bam file and index

#rm ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bam
#rm ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bam
#rm ${Mutect2_result}/SRR10351806_rg_added_sorted_dedupped_removed.realigned.bai
#rm ${Mutect2_result}/SRR10351807_rg_added_sorted_dedupped_removed.realigned.bai
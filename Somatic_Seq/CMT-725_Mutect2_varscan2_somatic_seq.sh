#!/bin/bash
#PBS -q batch
#PBS -N CMT-725-Mutect2-Varscan2-somatic_seq-filtering
#PBS -l nodes=1:ppn=1
#PBS -l walltime=30:00:00
#PBS -l mem=32gb


MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
reference='/work/szlab/Lab_shared_PanCancer/source'
Mutect2_result='/scratch/kh31516/Original_Mammary/store/Mutect2/CMT-725' 
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
Lofreq_folder='/scratch/kh31516/Original_Mammary/store/lofreq/CMT-725/loqfreq'
VarScan2_folder='/scratch/kh31516/Original_Mammary/store/varscan2/CMT-725/Vascan2'
somatic_seq_output_folder='/scratch/kh31516/Original_Mammary/store/somatic_seq/CMT-725'
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
script='/work/szlab/kh31516_Lab_Share_script'
#source_data='/scratch/kh31516/Original_Mammary/results/CMT-725'

: "
#mkdir Mutect_result
cd ${Mutect2_result}

ml SAMtools/1.9-GCC-8.2.0-2.31.1
Normal_sample=$(samtools view -H ${Mutect2_result}/SRR7780902_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${Mutect2_result}/SRR7780901_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)


## change the GATK version ##
module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8

gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${Mutect2_result}/SRR7780901_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${Mutect2_result}/SRR7780902_rg_added_sorted_dedupped_removed.realigned.bam \
--callable-depth 8 \
-normal $Normal_sample \
--panel-of-normals  ${MuTect2_source}/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${Mutect2_result}/CMT2-f1r2.tar.gz \
-L $GATK4/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${Mutect2_result}/CMT-725_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${Mutect2_result}/CMT2-f1r2.tar.gz -O ${Mutect2_result}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${Mutect2_result}/read-orientation-model.tar.gz \
-V ${Mutect2_result}/CMT-725_MuTect2_GATK4_noDBSNP.vcf \
-O ${Mutect2_result}/filtered-CMT-725_MuTect2_GATK4_noDBSNP.vcf


######### Varscan2 ###########
module load SAMtools/1.9-foss-2016b
ml VarScan/2.4.2-Java-1.8.0_144

cd ${VarScan2_folder}

#### Filter the somatic mutation derived from Varscan2 ####
time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter ${VarScan2_folder}/1mpileup-CMT-725_removed_realigned_Varscan.snp.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file ${VarScan2_folder}/CMT-725_Varscan2.vcf.snp_filtered.vcf

time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter ${VarScan2_folder}/1mpileup-CMT-725_removed_realigned_Varscan.indel.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file ${VarScan2_folder}/CMT-725_Varscan2.vcf.indel_filtered.vcf


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
--tumor-bam-file    ${Mutect2_result}/SRR7780901_rg_added_sorted_dedupped_removed.realigned.bam \
--normal-bam-file   ${Mutect2_result}/SRR7780902_rg_added_sorted_dedupped_removed.realigned.bam \
--mutect2-vcf       ${Mutect2_result}/filtered-CMT-725_MuTect2_GATK4_noDBSNP.vcf \
--lofreq-snv        ${Lofreq_folder}/CMT-725_lofreq.vcfsomatic_final.snvs.vcf.gz \
--lofreq-indel      ${Lofreq_folder}/CMT-725_lofreq.vcfsomatic_final.indels.vcf.gz \
--varscan-snv       ${VarScan2_folder}/CMT-725_Varscan2.vcf.snp_filtered.vcf \
--varscan-indel     ${VarScan2_folder}/CMT-725_Varscan2.vcf.indel_filtered.vcf 
"

####### filtering the consensus VCF files created by somaticseq with Dbsnp and Panel of Normals #########

cd $somatic_seq_output_folder

module load somaticseq/3.4.1_conda
module load BEDTools/2.28.0-foss-2018a
source activate ${SOMATICSEQROOT}


java -cp  ${script}/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${somatic_seq_output_folder}/Consensus.sSNV.vcf \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sSNV.vcf

java -cp  ${script}/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${somatic_seq_output_folder}/Consensus.sINDEL.vcf \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sINDEL.vcf


java -cp  ${script}/ DbSNP_filtering \
${MuTect2_source}/pon.vcf.txt \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sSNV.vcf \
${somatic_seq_output_folder}/PON_DbSNPFiltering_Consensus.sSNV.vcf

java -cp  ${script}/ DbSNP_filtering \
${MuTect2_source}/pon.vcf.txt \
${somatic_seq_output_folder}/DbSNPFiltering_Consensus.sINDEL.vcf \
${somatic_seq_output_folder}/PON_DbSNPFiltering_Consensus.sINDEL.vcf


### Change the filter status and prioritize by Mutect2 

cd $somatic_seq_output_folder

python ${script}/prioritizedByMutect2.py \
${somatic_seq_output_folder}/PON_DbSNPFiltering_Consensus.sSNV.vcf \
${somatic_seq_output_folder}/Mutect2_Prioritized_PON_DbSNPFiltering_Consensus.sSNV.vcf

python ${script}/prioritizedByMutect2.py \
${somatic_seq_output_folder}/PON_DbSNPFiltering_Consensus.sINDEL.vcf \
${somatic_seq_output_folder}/Mutect2_Prioritized_PON_DbSNPFiltering_Consensus.sINDEL.vcf
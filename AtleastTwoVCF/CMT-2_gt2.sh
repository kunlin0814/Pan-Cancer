#!/bin/bash
#PBS -q batch
#PBS -N CMT-2-Mutect2-Varscan2-strelka-somatic_seq-filtering
#PBS -l nodes=1:ppn=4
#PBS -l walltime=30:00:00
#PBS -l mem=32gb

bamFolder='/scratch/kh31516/Original_Mammary/results/CMT-2'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
reference='/work/szlab/Lab_shared_PanCancer/source'
Mutect2_result='/scratch/kh31516/Original_Mammary/store/Mutect2/CMT-2' 
GATK4='/work/szlab/Lab_shared_PanCancer/source/GATK4'
Lofreq_folder='/scratch/kh31516/Original_Mammary/store/lofreq/CMT-2/'
VarScan2_folder='/scratch/kh31516/Original_Mammary/store/varscan2/CMT-2/'
ConsensusFolder='/scratch/kh31516/Original_Mammary/store/ConsensusFolder/CMT-2'
script='/work/szlab/kh31516_Lab_Share_script'
strelkaOut='/scratch/kh31516/Original_Mammary/store/strelka_result/CMT-2-demo_somatic/results/variants'
#source_data='/scratch/kh31516/Original_Mammary/results/CMT-2'
Normal='SRR7780922'
Tumor='SRR7780923'

############################### Somatic Seq #################################
### Somatic Seq , needs the filtered Mutect2, filtered lofreq and filtered varscan2 vcf file, 
## The order of 3 VCF files cannot change, The orders are Mutect2 , lofreq, varscan2 (strelka)
## ( it will affect the further the analysis prioritized by Mutect2) ###

cd $bamFolder
module load picard/2.16.0-Java-1.8.0_144
## remove old index
rm ${bamFolder}/${Normal}_rg_added_sorted_dedupped_removed.realigned.bam.bai
rm ${bamFolder}/${Tumor}_rg_added_sorted_dedupped_removed.realigned.bam.bai
rm ${bamFolder}/${Normal}_rg_added_sorted_dedupped_removed.realigned.bai
rm ${bamFolder}/${Tumor}_rg_added_sorted_dedupped_removed.realigned.bai

## index tumor

java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex \
I=${bamFolder}/${Tumor}_rg_added_sorted_dedupped_removed.realigned.bam


## index normal
java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex \
I=${bamFolder}/${Normal}_rg_added_sorted_dedupped_removed.realigned.bam

mkdir -p $ConsensusFolder
cd  $ConsensusFolder
module load LoFreq/2.1.2-foss-2016b-Python-2.7.14
module load GATK/4.1.6.0-GCCcore-8.2.0-Java-1.8
ml SAMtools/1.9-GCC-8.2.0-2.31.1
ml VarScan/2.4.2-Java-1.8.0_144
module load somaticseq/3.4.1_conda
module load BEDTools/2.28.0-foss-2018a
source activate ${SOMATICSEQROOT}

## all the input vcf were filtered before somaticseq, 
# these input files didn't not go through dbsnp nor limit to CDS region

somaticseq_parallel.py \
--output-directory ${ConsensusFolder} \
--genome-reference ${reference}/canFam3.fa \
--inclusion-region ${reference}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_tab.bed \
--algorithm ada \
--threads 4 \
paired \
--tumor-bam-file    ${bamFolder}/${Tumor}_rg_added_sorted_dedupped_removed.realigned.bam \
--normal-bam-file   ${bamFolder}/${Normal}_rg_added_sorted_dedupped_removed.realigned.bam \
--mutect2-vcf       ${Mutect2_result}/filtered-CMT-2_MuTect2_GATK4_noDBSNP.vcf \
--lofreq-snv        ${Lofreq_folder}/CMT-2_lofreq.vcfsomatic_final.snvs.vcf.gz \
--lofreq-indel      ${Lofreq_folder}/CMT-2_lofreq.vcfsomatic_final.indels.vcf.gz \
--varscan-snv       ${VarScan2_folder}/CMT-2_Varscan2.vcf.snp_filtered.vcf \
--varscan-indel     ${VarScan2_folder}/CMT-2_Varscan2.vcf.indel_filtered.vcf \
--strelka-snv       ${strelkaOut}/somatic.snvs.vcf.gz \
--strelka-indel     ${strelkaOut}/somatic.indels.vcf.gz


#--strelka-snv       ${strelkaOut}/somatic.snvs.vcf.gz \
#--strelka-indel     ${strelkaOut}/somatic.indels.vcf.gz


####### filtering the consensus VCF files created by somaticseq with Dbsnp #########


java -cp  ${script}/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${ConsensusFolder}/Consensus.sSNV.vcf \
${ConsensusFolder}/CMT-2-DbSNPFiltering_Consensus.sSNV.vcf

java -cp  ${script}/ DbSNP_filtering \
$reference/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${ConsensusFolder}/Consensus.sINDEL.vcf \
${ConsensusFolder}/CMT-2-DbSNPFiltering_Consensus.sINDEL.vcf



### Grep At least two software ###

cd $ConsensusFolder

python ${script}/grepAtLeastTwoVcf.py \
${ConsensusFolder}/CMT-2-DbSNPFiltering_Consensus.sSNV.vcf \
${ConsensusFolder}/AtLeastTwo_CMT-2-DbSNPFiltering_Consensus.sSNV.vcf


python ${script}/grepAtLeastTwoVcf.py \
${ConsensusFolder}/CMT-2-DbSNPFiltering_Consensus.sINDEL.vcf \
${ConsensusFolder}/AtLeastTwo_CMT-2-DbSNPFiltering_Consensus.sINDEL.vcf

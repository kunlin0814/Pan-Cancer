#!/bin/bash
#PBS -N VarScan2-lofreq-CMT-2
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -q batch


cd $PBS_O_WORKDIR

result='/scratch/kh31516/Original_Mammary/results/CMT-2' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
VarScan_outfolder='/scratch/kh31516/Original_Mammary/results/CMT-2/Vascan2'
lofreq_outfolder='/scratch/kh31516/Original_Mammary/results/CMT-2/loqfreq'
Lofeq_source='/work/szlab/Lab_shared_PanCancer/source/Lofreq'
CDS_script='/work/szlab/kh31516_Lab_Share_script'

#mkdir $VarScan_outfolder
cd $result


module load SAMtools/1.9-foss-2016b
ml VarScan/2.4.2-Java-1.8.0_144
#### Sort Bam File ####

# Normal
#samtools sort $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam -o $result/sort-SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam

# Tumor
#samtools sort $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam -o $result/sort-SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam


## Create mpileup files Normal and tumor sample together

#samtools mpileup -f $reference/canFam3.fa -q 1 -B $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
#$result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam > $VarScan_outfolder/SRR7780922-SRR7780923.mpileup

#samtools mpileup -f $reference/canFam3.fa -q 1 -B $result/sort-SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
#$result/sort-SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam > $VarScan_outfolder/SRR7780922-SRR7780923.mpileup


## Their parameters
#--min-coverage 8 
#--min-coverage-normal 6 
#--min-coverage-tumor 8 
#--min-reads2 2 
#--min-avg-qual 15 
#--min-var-freq 0.08 
#--min-freq-for-hom 0.75 
#--tumor-purity 1.0 
#--strand-filter 1 
#--somatic-p-value 0.05 
#--output-vcf 1

cd $VarScan_outfolder

##### VarScan2 #####

: "
time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somatic \
$VarScan_outfolder/SRR7780922-SRR7780923.mpileup \
$VarScan_outfolder/1mpileup-CMT-2_removed_realigned_Varscan \
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
"

#### Filter the somatic mutation derived from Varscan2 ####
time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter $VarScan_outfolder/1mpileup-CMT-2_removed_realigned_Varscan.snp.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file $VarScan_outfolder/1mpileup-CMT-2_Varscan.vcf.snp_filtered.vcf

time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar somaticFilter $VarScan_outfolder/1mpileup-CMT-2_removed_realigned_Varscan.indel.vcf \
--min-coverage 8 \
-–min-reads2 2 \
-–min-var-freq 0.08 \
--p-value 0.05 \
-–min-avg-qual 15 \
--output-file $VarScan_outfolder/1mpileup-CMT-2_Varscan.vcf.indel_filtered.vcf


#### This command limits variants in a file to a set of positions or regions ####

## need to pass first and then limits

cat $VarScan_source_folder/1mpileup-CMT-2_Varscan.vcf.snp_filtered.vcf | awk '{if ($7=="PASS") {print $0}}' > $VarScan2_folder/1mpileup-CMT-2_Varscan.vcf.snp_filtered.vcf-PASS.vcf

cat $VarScan_source_folder/HSA_1_Varscan2.vcf.indel_filtered.vcf | awk '{if ($7=="PASS") {print $0}}' > $VarScan2_folder/HSA_1_Varscan2.vcf.indel_filtered.vcf-PASS.vcf




time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar limit $VarScan_outfolder/1mpileup-CMT-2_removed_realigned_Varscan.snp.vcf \
--regions-file $reference/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_tab \
--output-file $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.snp.vcf


time java -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar limit $VarScan_outfolder/1mpileup-CMT-2_Varscan.vcf.indel_filtered.vcf-PASS.vcf \
--regions-file $reference/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_tab \
--output-file $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.indel.vcf

### grep PASS ###

cat $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.snp.vcf | grep -w "PASS" >  $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.snp.vcf.PASS

cat $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.indel.vcf |grep -w "PASS" > $VarScan_outfolder/CDS-1mpileup-CMT-2_removed_realigned_Varscan.indel.vcf.PASS


rm $VarScan_outfolder/CDS-1mpileup-CMT-2_Varscan.vcf.indel_filtered.vcf-PASS.vcf
rm $VarScan_outfolder/CDS-1mpileup-CMT-2_Varscan.vcf.snp_filtered.vcf-PASS.vcf

: "
###### LoFreq ######
mkdir $lofreq_outfolder
cd $lofreq_outfolder

module load LoFreq/2.1.2-foss-2016b-Python-2.7.14

lofreq somatic -n $result/SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam \
-t $result/SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam \
-f $reference/canFam3.fa \
--threads 4 -o $lofreq_outfolder/CMT-2_lofreq.vcf \
-d $Lofeq_source/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf.gz \
--call-indels \
--min-cov 7 \
--verbose

### Limit the variants in the CDS region ###
python3 $CDS_script/lofreq_limitCDS.py \
$lofreq_outfolder/CMT-2_lofreq.vcfsomatic_final_minus-dbsnp.snvs.vcf.gz \
$reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
$lofreq_outfolder/CMT-2_SNV

### Limit the indels in the CDS region ###
python3 $CDS_script/lofreq_limitCDS.py \
$lofreq_outfolder/CMT-2_lofreq.vcfsomatic_final_minus-dbsnp.indels.vcf.gz \
$reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
$lofreq_outfolder/CMT-2_indels

#rm $result/SRR7780922_removed.realigned.mpileup
#rm $result/SRR7780923_removed.realigned.mpileup
rm $VarScan_outfolder/SRR7780922-SRR7780923.mpileup
rm $VarScan_outfolder/1mpileup-CMT-2_Varscan.vcf.snp_filtered.vcf
#rm $result/sort-SRR7780922_rg_added_sorted_dedupped_removed.realigned.bam
#rm $result/sort-SRR7780923_rg_added_sorted_dedupped_removed.realigned.bam
"
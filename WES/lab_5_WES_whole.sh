#!/bin/bash
#SBATCH --job-name=lab5_WES         # Job name (lab5_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=40G                   # Job memory limit (10 GB)
#SBATCH --time=40:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=lab5_WES.%j.out    # Standard output log
#SBATCH --error=lab5_WES.%j.err     # Standard error log


################################################################################################################################
###### complete pipeline (picard, GATK, annova, Mutect, Mutect2,Strelka, Germline mutation, Depth of coverage, Callablebases) ##########
##### WGS analysis needs to change the interval region #####


thread=4
Bioproject='Lab_WES_MC'
Normal_Run='5M-SL7497'
Tumor_Run='5T-SL7496'
SampleName='lab_5'
results='/scratch/yf94402/Pan_cancer/data/' #
relaign_output='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/relaign_output/'${SampleName}'' #
reference='/work/szlab/Lab_shared_PanCancer/source' 
dataN='/scratch/kh31516/Original_Mammary/data/'${Bioproject}'/'${SampleName}'/'${Normal_Run}
dataT='/scratch/kh31516/Original_Mammary/data/'${Bioproject}'/'${SampleName}'/'${Tumor_Run}
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
script='/home/kh31516/kh31516_Lab_Share_script'
MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
MuTect_out='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/Mutect/'${SampleName}''
MuTect2_out='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/Mutect2/'${SampleName}''
Germline_out='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/Germline/'${SampleName}''
DepthOfCoverage='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/DepthOfCoverage/'${SampleName}''
strelka_out='/scratch/kh31516/Original_Mammary/store/'${Bioproject}'/strelka/'${SampleName}''
###sraRunTable='/scratch/kh31516/Pan_cancer/PRJNA247493/source/PRJNA247493_SraRunTable.txt'# path



mkdir -p ${MuTect_out}
mkdir -p ${MuTect2_out}
mkdir -p ${Germline_out}
mkdir -p ${DepthOfCoverage}
mkdir -p $strelka_out
mkdir -p $relaign_output
####### GATK Realign #######
# Generating interval file for sort.bam
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
ml Anaconda3/2020.02

cd ${MuTect_out}
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar $EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# add gene name
python2 $script/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



### Sanger 5 steps filtering ###
# 5 Steps filtering
grep -w KEEP ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${MuTect_out}/${SampleName}_PASS.stat

python $script/Filter_MutectStat_5steps.py \
${MuTect_out}/${SampleName}_PASS.stat \
${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
${MuTect_out}/${SampleName}_vaf_before.txt \
${MuTect_out}/${SampleName}_vaf_after.txt \
${MuTect_out}/${SampleName}_whyout.txt 

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python2 $script/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


############################## Germline mutation preparation Start ######################################
# GATK 3 
ml GATK/3.8-1-Java-1.8.0_144

cd $Germline_out

# Variant calling
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/canFam3.fa -I \
${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam -dontUseSoftClippedBases \
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
-I ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Normal_Run}_table.txt \
-o ${Germline_out}/${Normal_Run}_callable_status.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T CallableLoci \
-R ${reference}/canFam3.fa \
-I ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-summary ${Germline_out}/${Tumor_Run}_table.txt \
-o ${Germline_out}/${Tumor_Run}_callable_status.bed 

############################## Germline mutation preparation End ######################################


#### GATK DepthofCoverage ####
cd ${DepthOfCoverage}
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Normal_Run}_DepthofCoverage_CDS.bed 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T DepthOfCoverage \
-R ${reference}/canFam3.fa \
-I ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--minBaseQuality 10 \
--minMappingQuality 10 \
-L ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
-o ${DepthOfCoverage}/${Tumor_Run}_DepthofCoverage_CDS.bed 


####### Mutect2 #######

cd ${MuTect2_out}
module load SAMtools/1.9-GCC-8.3.0
module load GATK/4.1.6.0-GCCcore-8.3.0-Java-1.8

Normal_sample=$(samtools view -H ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)
Tumor_sample=$(samtools view -H ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam | grep '@RG' |awk -F "SM:" '{print $2}' | cut -f1)

## Mutect2 follows the Glioma PanCancer parameters
gatk Mutect2 \
-R ${reference}/canFam3.fa \
-I ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-I ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
-normal $Normal_sample \
--callable-depth 8 \
--panel-of-normals  $MuTect2_source/pon.vcf.gz \
--initial-tumor-lod 2.0 --normal-lod 2.2 --tumor-lod-to-emit 3.0 --pcr-indel-model CONSERVATIVE \
--f1r2-tar-gz ${MuTect2_out}/${SampleName}-f1r2.tar.gz \
--dontUseSoftClippedBases true \
-L ${MuTect2_source}/Uniq-Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval.list \
-O ${MuTect2_out}/${SampleName}_MuTect2_GATK4_noDBSNP.vcf


## pass this raw data to LearnReadOrientationModel:
gatk LearnReadOrientationModel -I ${MuTect2_out}/${SampleName}-f1r2.tar.gz -O ${MuTect2_out}/read-orientation-model.tar.gz

# filtered Mutect2 files

gatk FilterMutectCalls -R ${reference}/canFam3.fa \
-ob-priors ${MuTect2_out}/read-orientation-model.tar.gz \
-V ${MuTect2_out}/${SampleName}_MuTect2_GATK4_noDBSNP.vcf \
-O ${MuTect2_out}/filtered-${SampleName}_MuTect2_GATK4_noDBSNP.vcf

cd ${MuTect2_out}

### filter Mutect2 #####

## with DbSNP
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/filtered-${SampleName}_MuTect2_GATK4_noDBSNP.vcf \
${MuTect2_out}/DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf


## with PON
java -cp  $script/ DbSNP_filtering \
${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
${MuTect2_out}/DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf


### Mutect2 5 steps filtering

cd ${MuTect2_out}

ml Anaconda3/2020.02
#source activate py35

#awk '$7 == "PASS" {print $0}' ${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf > ${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf-PASS

# 5 Steps filtering
python $script/Mutect2_5Steps_filtering.py \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/PON_DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf \
${MuTect2_out}/${SampleName}_VAF_Before.txt \
${MuTect2_out}/${SampleName}_VAF_After.txt


## Strelka
cd ${results}
module load SAMtools/1.9-GCC-8.3.0
#samtools index ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam
#samtools index ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam

${script}/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${relaign_output}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--tumorBam ${relaign_output}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--referenceFasta ${reference}/canFam3.fa \
--runDir $strelka_out

# your result will be generated in demo_somatic/results/variants
$strelka_out/runWorkflow.py -m local -j 20


## limit strelka result into CDS region

python2 $script/Limit_vcf_to_CDS.py $strelka_out/results/variants/somatic.indels.vcf.gz $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list


### Annovar for strelka
cd $strelka_out/results/

####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS > $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS

# annovar input preparation
perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS > \
$strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput

# annovar annotate
perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput \
$reference/annovar_CanFam3.1.99.gtf

# add gene name
python2 $script/Add_GeneName_N_Signature.py $strelka_out/results/variants/somatic.indels.vcf_canFam3.1.99_CDS-PASS-avinput.exonic_variant_function \
$reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


: "
####### remove unneeded files #######
rm ${results}/*_sorted_dedupped_removed.bam
rm ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm ${results}/*.bai
rm ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS ${Germline_out}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf-PASS-avinput
rm ${results}/${Normal_Run}.sam ${results}/${Tumor_Run}.sam 
rm ${results}/header_N ${results}/header_T
rm ${results}/${Normal_Run}_rg_added_sorted.bam ${results}/${Tumor_Run}_rg_added_sorted.bam
rm ${results}/${Normal_Run}-cleaned.sam ${results}/${Tumor_Run}-cleaned.sam
rm $dataN/*.fastq.gz $dataT/*.fastq.gz
rm ${results}/*.sai
rm ${MuTect2_out}/DbSNP_filtered-${SampleName}_MuTect2_GATK4.vcf
rm ${MuTect2_out}/read-orientation-model.tar.gz
rm ${MuTect2_out}/${SampleName}-f1r2.tar.gz
"

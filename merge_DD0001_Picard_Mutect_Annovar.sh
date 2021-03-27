#!/bin/bash
#SBATCH --job-name=DD0001_merge_Mutect_WES         # Job name (DD0001_WES)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task (4)
#SBATCH --mem=60G                   # Job memory limit (10 GB)
#SBATCH --time=100:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=DD0001_merge_Mutect_WES.%j.out    # Standard output log
#SBATCH --error=DD0001_merge_Mutect_WES.%j.err     # Standard error log

thread=4
Bioproject='PRJNA552034'
Normal_Run1='ERR1681499'
Normal_Run2='ERR1681454'
Tumor_Run1='ERR1681500'
Tumor_Run2='ERR1681455'
final_Normal_Run='ERR1681499-ERR1681454'
final_Tumor_Run='ERR1681500-ERR1681455'
SampleName='DD0001a'
results='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/results/'${SampleName} #
reference='/work/szlab/Lab_shared_PanCancer/source' 
MuTect_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/'${SampleName}''
# MuTect2_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect2/'${SampleName}''
# MuTect2_source='/work/szlab/Lab_shared_PanCancer/source/Mutect2'
# Germline_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/Germline/'${SampleName}''
# DepthOfCoverage='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/DepthOfCoverage/'${SampleName}''
strelka_out='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/Re-run_OM_Picard_mutect/Mutect/strelka/'${SampleName}''
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'
scripts='/home/kh31516/kh31516_Lab_Share_script'


module load SAMtools/1.9-GCC-8.3.0
module load SRA-Toolkit/2.9.6-1-centos_linux64
module load picard/2.21.6-Java-11
module load Perl/5.26.1-GCCcore-6.4.0
# ml BWA/0.7.17-GCC-8.3.0

mkdir -p $MuTect_out

# ####### Convert sam to bam file #######

# cd ${results}
# module load SAMtools/1.9-GCC-8.3.0

# samtools view -bS@ 4 ${results}/${Normal_Run1}.sam > ${results}/${Normal_Run1}.bam
# samtools view -bS@ 4 ${results}/${Normal_Run2}.sam > ${results}/${Normal_Run2}.bam
# samtools view -bS@ 4 ${results}/${Tumor_Run1}.sam > ${results}/${Tumor_Run1}.bam
# samtools view -bS@ 4 ${results}/${Tumor_Run2}.sam > ${results}/${Tumor_Run2}.bam

# ## merge_bam ##
cd ${results}
samtools merge ${results}/${final_Normal_Run}.bam ${results}/${Normal_Run1}.bam ${results}/${Normal_Run2}.bam
samtools merge ${results}/${final_Tumor_Run}.bam ${results}/${Tumor_Run1}.bam ${results}/${Tumor_Run2}.bam 

####### Picard #######
# get header information
samtools view -H ${results}/${final_Normal_Run}.bam > $results/header_N
samtools view -H ${results}/${final_Tumor_Run}.bam > $results/header_T

# exclude unmapped reads based on FLAG
samtools view $results/${final_Normal_Run}.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $results/${final_Normal_Run}-cleaned.sam
cat $results/header_N $results/${final_Normal_Run}-cleaned.sam > $results/fooN
mv $results/fooN $results/${final_Normal_Run}-cleaned.sam
samtools view $results/${final_Tumor_Run}.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $results/${final_Tumor_Run}-cleaned.sam
cat $results/header_T $results/${final_Tumor_Run}-cleaned.sam > $results/fooT
mv $results/fooT $results/${final_Tumor_Run}-cleaned.sam

####### Picard #######
# picard sort
cd ${results}
module load picard/2.21.6-Java-11

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I=${results}/${final_Normal_Run}-cleaned.sam \
-O=${results}/${final_Normal_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${final_Normal_Run}_normal

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I=${results}/${final_Tumor_Run}-cleaned.sam \
-O=${results}/${final_Tumor_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${final_Tumor_Run}_tumor

# picard MarkDuplicates
java -jar $EBROOTPICARD/picard.jar  MarkDuplicates -I=${results}/${final_Normal_Run}_rg_added_sorted.bam \
-O=${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${final_Normal_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar  MarkDuplicates -I=${results}/${final_Tumor_Run}_rg_added_sorted.bam \
-O=${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${final_Tumor_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

# ####### Remove unneeded files #######
# rm $results/${final_Normal_Run}.sam $results/ERR1681500.sam 
# rm $results/header_N $results/header_T
# rm $results/${final_Normal_Run}_rg_added_sorted.bam $results/${final_Tumor_Run}_rg_added_sorted.bam
# rm $results/${final_Normal_Run}-cleaned.sam $results/${final_Tumor_Run}-cleaned.sam

ml GATK/3.8-1-Java-1.8.0_144
####### GATK Realign #######
# Generating interval file for sort.bam

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt $thread \
-o ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar $EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect \ \
--reference_sequence ${reference}/canFam3.fa \
--dbsnp ${reference}/DbSNP_canFam3_version151-DogSD_Feb2020_V4.vcf \
--defaultBaseQualities 30 \
--intervals ${reference}/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal ${results}/${final_Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor ${results}/${final_Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam \
--out ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $results/ERR1681500_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $results/${final_Normal_Run}_rg_added_sorted_dedupped_removed.realigned \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf



####### Annovar #######
# Extract PASS records from vcf
module load Perl/5.26.1-GCCcore-6.4.0

awk '$7 == "PASS" {print $0}' ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index

# add gene name
python2 $scripts/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt



### Sanger 5 steps filtering ###
# 5 Steps filtering

## With 5 steps filtering 

ml Anaconda3/2020.02
grep -w KEEP ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > ${MuTect_out}/${SampleName}_PASS.stat

python $scripts/Filter_MutectStat_5steps.py \
$mutect_folder/${SampleName}/${SampleName}_PASS.stat \
$mutect_folder/${SampleName}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
$mutect_folder/${SampleName}/${SampleName}_vaf_before.txt \
$mutect_folder/${SampleName}/${SampleName}_vaf_after.txt \
$mutect_folder/${SampleName}/${SampleName}_whyout.txt \
$SampleName 

## 5 steps filtering created XXX_filterMut file 

# annovar input preparation
perl $annovar_index/convert2annovar.pl -format vcf4old ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $annovar_index/annotate_variation.pl --buildver canFam3 ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $annovar_index

# add gene name
python2 $scripts/Add_GeneName_N_Signature.py ${MuTect_out}/${SampleName}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt


# ####### remove unneeded files #######
# rm $results/*_sorted_dedupped_removed.bam
# rm $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $results/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
# rm $results/*.bai
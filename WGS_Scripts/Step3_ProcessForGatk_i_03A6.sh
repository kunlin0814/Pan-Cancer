#!/bin/bash
#SBATCH --job-name=i_03A6_processGATK         # Job name (i_03A6)
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

####### Convert sam to bam file #######

cd ${results}
ml SAMtools/1.9-GCC-8.3.0

samtools view -bS ${results}/${Tumor_Run}.sam > ${results}/${Tumor_Run}.bam
samtools view -bS ${results}/${Normal_Run}.sam > ${results}/${Normal_Run}.bam

# get header information
#grep '@SQ\|@PG' ${results}/${Normal_Run}.sam > ${results}/header_N
#grep '@SQ\|@PG' ${results}/${Tumor_Run}.sam > ${results}/header_T
samtools view -H ${results}/${Normal_Run}.bam > ${results}/header_N
samtools view -H ${results}/${Tumor_Run}.bam > ${results}/header_T

# exclude unmapped reads based on FLAG
cat ${results}/${Normal_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Normal_Run}-cleaned.sam
cat ${results}/header_N ${results}/${Normal_Run}-cleaned.sam > ${results}/fooN
mv ${results}/fooN ${results}/${Normal_Run}-cleaned.sam

cat ${results}/${Tumor_Run}.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > ${results}/${Tumor_Run}-cleaned.sam
cat ${results}/header_T ${results}/${Tumor_Run}-cleaned.sam > ${results}/fooT
mv ${results}/fooT ${results}/${Tumor_Run}-cleaned.sam


####### Picard #######
# picard sort
cd ${results}
module load picard/2.21.6-Java-11

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Normal_Run}-cleaned.sam \
O=${results}/${Normal_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Normal_Run}_normal

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${results}/${Tumor_Run}-cleaned.sam \
O=${results}/${Tumor_Run}_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${Tumor_Run}_tumor

# picard MarkDuplicates
java -jar $EBROOTPICARD/picard.jar  MarkDuplicates I=${results}/${Normal_Run}_rg_added_sorted.bam \
O=${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Normal_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true

java -jar $EBROOTPICARD/picard.jar  MarkDuplicates I=${results}/${Tumor_Run}_rg_added_sorted.bam \
O=${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT M=${results}/${Tumor_Run}-output.metrics REMOVE_SEQUENCING_DUPLICATES=true


####### GATK Realign #######
# Generating interval file for sort.bam
ml GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-nt 4 \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

#### IndelRealigner realign ######

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Normal_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/canFam3.fa \
-I ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam \
-targetIntervals ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.bam.intervals \
-o ${results}/${Tumor_Run}_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores

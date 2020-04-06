#!/bin/bash
#SBATCH --job-name=GATK_NCT_MuTect # Job name
#SBATCH --partition=batch # Partition (queue) name
#SBATCH --ntasks=1 # Single task job
#SBATCH --cpus-per-task=8 # Number of cores per task
#SBATCH --mem=35gb # Total memory for job
#SBATCH --time=150:00:00 # Time limit hrs:min:sec
#SBATCH --output=RestOf.%j # Standard output and error log
#SBATCH --mail-user=kh31516@uga.edu
#SBATCH --mail-type=END

###### complete pipeline (including callable base) ##########
cd $SLURM_SUBMIT_DIR

data_source='/home/kh31516/8211/HW6/results-3' #
result='/home/kh31516/8211/HW6/results-4'
reference='/home/kh31516/8211/HW6/source' 
data='/home/kh31516/8211/HW6/data'

cd ${result}

ml GATK/3.8-0-Java-1.8.0_144

### Generating recalibration table ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${data_source}/SRR10993129_rg_added_sorted_dedupped_removed.bam \
-nct 8 \
-knownSites ${reference}/dbsnp_146.hg38.vcf -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed.bam.recal.table

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${data_source}/SRR10993128_rg_added_sorted_dedupped_removed.bam \
-nct 8 \
-knownSites ${reference}/dbsnp_146.hg38.vcf -o ${result}/SRR10993128_rg_added_sorted_dedupped_removed.bam.recal.table

### Calibrating data ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${data_source}/SRR10993129_rg_added_sorted_dedupped_removed.bam \
-nct 8 \
-BQSR ${result}/SRR10993129_rg_added_sorted_dedupped_removed.bam.recal.table -o ${result}/SRR10993128_rg_added_sorted_dedupped_removed_recal.bam

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T PrintReads -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${data_source}/SRR10993128_rg_added_sorted_dedupped_removed.bam \
-nct 8 \
-BQSR ${result}/SRR10993128_rg_added_sorted_dedupped_removed.bam.recal.table -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam

### RealignerTargetCreator ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam.intervals \
-nt 8 \
--allow_potentially_misencoded_quality_scores

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam.intervals \
-nt 8 \
--allow_potentially_misencoded_quality_scores

### Indel realigning ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam -targetIntervals ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam.intervals \
-o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_realigned.bam --allow_potentially_misencoded_quality_scores

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam -targetIntervals ${result}/SRR10993129_rg_added_sorted_dedupped_removed_recal.bam.intervals \
-o ${result}/SRR10993128_rg_added_sorted_dedupped_removed_realigned.bam --allow_potentially_misencoded_quality_scores


#### HaplotypeCaller ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993128_rg_added_sorted_dedupped_removed_realigned.bam -o ${result}/SRR10993128_rg_added_sorted_dedupped_removed_realigned.g.vcf \
-ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/Homo_sapiens_assembly38.fasta \
-I ${result}/SRR10993129_rg_added_sorted_dedupped_removed_realigned.bam -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_realigned.g.vcf \
-ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000

### Variant filtering ###
time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/Homo_sapiens_assembly38.fasta \
-V ${result}/SRR10993128_rg_added_sorted_dedupped_removed_realigned.g.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" \
-o ${result}/SRR10993128_rg_added_sorted_dedupped_removed_realigned_filtered.g.vcf

time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/Homo_sapiens_assembly38.fasta \
-V ${result}/SRR10993129_rg_added_sorted_dedupped_removed_realigned.g.vcf \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed_realigned_filtered.g.vcf



# lab version Variant calling
#time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/Homo_sapiens_assembly38.fasta -I ${result}/SRR10993129_rg_added_sorted_dedupped_removed.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${result}/SRR10993128_rg_added_sorted_dedupped_removed.realigned.bam.vcf
#time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference}/Homo_sapiens_assembly38.fasta -I ${result}/SRR10993128_rg_added_sorted_dedupped_removed.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed.realigned.bam.vcf

# Variant filtering
#time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/Homo_sapiens_assembly38.fasta -V ${result}/SRR10993128_rg_added_sorted_dedupped_removed.realigned.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${result}/SRR10993128_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf
#time java -Xmx4g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantFiltration -R ${reference}/Homo_sapiens_assembly38.fasta -V ${result}/SRR10993129_rg_added_sorted_dedupped_removed.realigned.bam.vcf -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${result}/SRR10993129_rg_added_sorted_dedupped_removed.realigned.bam.filter.vcf

## Mutect

module load MuTect/1.1.7-Java-1.7.0_80

java -Xmx4g -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence /home/kh31516/8211/HW6/source/Homo_sapiens_assembly38.fasta \
--dbsnp /home/kh31516/8211/HW6/source/dbsnp_146.hg38.vcf \
--intervals /home/kh31516/8211/HW6/source/NISTIntegratedCalls.hg38.interval_list \
--input_file:normal /home/kh31516/8211/HW6/results-3/SRR10993128_rg_added_sorted_dedupped_removed.bam \
--input_file:tumor /home/kh31516/8211/HW6/results-3/SRR10993129_rg_added_sorted_dedupped_removed.bam \
--out /home/kh31516/8211/HW6/results-4/Human_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file /home/kh31516/8211/HW6/results-4/GATK-Human_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf /home/kh31516/8211/HW6/results-4/GATK-Human_rg_added_sorted_dedupped_removed.MuTect.vcf


module load MuTect/1.1.7-Java-1.7.0_80

java -Xmx4g -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence ${reference}/Homo_sapiens_assembly38.fasta \
--dbsnp ${reference}/dbsnp_146.hg38.vcf \
--intervals ${reference}/CDS-Homo_sapiens.GRCh38.99.chr.interval_list \
--input_file:normal ${result}/SRR10993128_rg_added_sorted_dedupped_removed.bam \
--input_file:tumor ${result}/SRR10993129_rg_added_sorted_dedupped_removed.bam \
--out ${result}/Our_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file ${result}/Human_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf ${result}/Our_rg_added_sorted_dedupped_removed.MuTect.vcf

#### Performing joint analysis ###
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference/hg38.fa -V $result/SRR10993128_rg_added_sorted_dedupped_removed_realigned_recal.g.vcf -V $result/SRR10993129_rg_added_sorted_dedupped_removed_realigned_recal.g.vcf -o $result/SRR10993128_SRR10993129_joint.vcf

#### Building Gaussian mixture model ###
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $reference/hg38.fa -input $result/SRR10993128_SRR10993129_joint.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
#-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz -an DP -an QD -an MQRankSum -an FS -mode SNP -recalFile $result/SRR10993128_SRR10993129_joint.recal \
#-tranchesFile $result/SRR10993128_SRR10993129_joint.tranches -rscriptFile $result/recal.plots.R

#### Applying recalibration ###
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R $reference/hg38.fa -input $result/SRR10993128_SRR10993129_joint.vcf -mode SNP -recalFile $result/SRR10993128_SRR10993129_joint.recal -tranchesFile $result/SRR10993128_SRR10993129_joint.tranches -o $result/SRR10993128_SRR10993129_joint_recal.vcf \
#-ts_filter_level 99.0
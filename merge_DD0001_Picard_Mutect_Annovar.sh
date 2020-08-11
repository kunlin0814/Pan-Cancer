#PBS -S /bin/bash
#PBS -q batch
#PBS -N DD0001-merge
#PBS -l nodes=1:ppn=4
#PBS -l walltime=100:00:00
#PBS -l mem=60gb
#PBS -M kh31516@uga.edu 
#PBS -m ae
cd /scratch/kh31516/Melanoma_merge/results/
result='/scratch/kh31516/Melanoma_merge/results/DD0001' #
reference='/scratch/kh31516/Melanoma_source' 


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
module load Perl/5.26.1-GCCcore-6.4.0

mkdir $result

## merge_bam ##
samtools merge ${result}/ERR1681499-ERR1681454.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-1/ERR1681499.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-2/ERR1681454.bam
samtools merge ${result}/ERR1681500-ERR1681455.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-1/ERR1681500.bam /scratch/kh31516/Melanoma_new/WES/results/DD0001-2/ERR1681455.bam  


####### Picard #######
# get header information
samtools view -H ${result}/ERR1681499-ERR1681454.bam > $result/header_N
samtools view -H ${result}/ERR1681500-ERR1681455.bam > $result/header_T

# exclude unmapped reads based on FLAG
samtools view $result/ERR1681499-ERR1681454.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/ERR1681499-ERR1681454-cleaned.sam
cat $result/header_N $result/ERR1681499-ERR1681454-cleaned.sam > $result/fooN
mv $result/fooN $result/ERR1681499-ERR1681454-cleaned.sam
samtools view $result/ERR1681500-ERR1681455.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/ERR1681500-ERR1681455-cleaned.sam
cat $result/header_T $result/ERR1681500-ERR1681455-cleaned.sam > $result/fooT
mv $result/fooT $result/ERR1681500-ERR1681455-cleaned.sam

# picard sort
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/ERR1681499-ERR1681454-cleaned.sam O=$result/ERR1681499-ERR1681454_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/ERR1681500-ERR1681455-cleaned.sam O=$result/ERR1681500-ERR1681455_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/ERR1681499-ERR1681454_rg_added_sorted.bam O=$result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/ERR1681499-output.metrics REMOVE_SEQUENCING_DUPLICATES=true
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/ERR1681500-ERR1681455_rg_added_sorted.bam O=$result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/ERR1681500-output.metrics REMOVE_SEQUENCING_DUPLICATES=true





####### Remove unneeded files #######
rm $result/ERR1681499-ERR1681454.sam $result/ERR1681500.sam 
rm $result/header_N $result/header_T
rm $result/ERR1681499-ERR1681454_rg_added_sorted.bam $result/ERR1681500-ERR1681455_rg_added_sorted.bam
rm $result/ERR1681499-ERR1681454-cleaned.sam $result/ERR1681500-ERR1681455-cleaned.sam





####### GATK Realign #######
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.bam -o $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.bam -o $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

# realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.bam.intervals -o $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.bam.intervals -o $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores





######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/ERR1681500-ERR1681455_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/DD0001_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/DD0001_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/ERR1681500_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/ERR1681499-ERR1681454_rg_added_sorted_dedupped_removed.realigned \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf





####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference






####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
rm $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/DD0001_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm $result/*.bai
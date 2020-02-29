#PBS -S /bin/bash
#PBS -N ND10_441
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=50gb
#PBS -q batch
#PBS -M kh31516@uga.edu 
#PBS -m ae

result='/scratch/kh31516/Melanoma/WGS/results/ND10_441' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
dataN='/scratch/kh31516/Melanoma/WGS/data/SAMN07376273'
dataT='/scratch/kh31516/Melanoma/WGS/data/SAMN07376267'

#module load MuTect/1.1.7-Java-1.7.0_80
module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0




cd $result

samtools index $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned.bam
samtools index $result/SAMN07376267_rg_added_sorted_dedupped_removed.realigned.bam


######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
#module load MuTect/1.1.7-Java-1.7.0_80

#time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
#--reference_sequence $reference/canFam3.fa \
#--dbsnp $reference/dbsnp_9615.vcf \
#--intervals $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#--input_file:normal $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned.bam \
#--input_file:tumor $result/SAMN07376267_rg_added_sorted_dedupped_removed.realigned.bam \
#--out $result/ND10_441_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
#--coverage_file $result/ND10_441_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
#--vcf $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf

########  Mutect cover whole region #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.interval_list \
--input_file:normal $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/SAMN07376267_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf


####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/ERR1681500_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf



####### Annovar #######
## Annovar to cover new region ##
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/whole_region_ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference






####### remove unneeded files #######
#rm $result/*_sorted_dedupped_removed.bam
#rm $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
#rm $result/*.bai
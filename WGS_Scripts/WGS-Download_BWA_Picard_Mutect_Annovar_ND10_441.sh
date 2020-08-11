#PBS -S /bin/bash
#PBS -N ND10_441
#PBS -l walltime=660:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=60gb
#PBS -q batch
#PBS -M kh31516@uga.edu 
#PBS -m ae

result='/scratch/kh31516/Melanoma/WGS/results/ND10_441' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
dataN='/scratch/kh31516/Melanoma/WGS/data/SAMN07376273'
dataT='/scratch/kh31516/Melanoma/WGS/data/SAMN07376267'


module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0


####### Download normal #######
## normal
mkdir $result
mkdir $dataN
cd $dataN
fasterq-dump SRR5859760 -t /scratch/kh31516/Gtex/tmp/ -e 6 --split-files
fasterq-dump SRR5867236 -t /scratch/kh31516/Gtex/tmp/ -e 6 --split-files

gzip $dataN/SRR5859760_1.fastq
gzip $dataN/SRR5859760_2.fastq

gzip $dataN/SRR5867236_1.fastq
gzip $dataN/SRR5867236_2.fastq


#fastq-dump --split-files --gzip ERR1681456

####### Download tumor #######
mkdir $dataT
cd $dataT
fasterq-dump SRR5859757 -t /scratch/kh31516/Gtex/tmp/ -e 6 --split-files
fasterq-dump SRR5867229 -t /scratch/kh31516/Gtex/tmp/ -e 6 --split-files

gzip $dataT/SRR5859757_1.fastq
gzip $dataT/SRR5859757_2.fastq

gzip $dataT/SRR5867229_1.fastq
gzip $dataT/SRR5867229_2.fastq

####### BWA Mapping #######
cd $result
# Tumor
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR5859757_1.fastq.gz > $result/SRR5859757_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR5859757_2.fastq.gz > $result/SRR5859757_2.sai

time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR5867229_1.fastq.gz > $result/SRR5867229_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataT/SRR5867229_2.fastq.gz > $result/SRR5867229_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR5859757_1.sai \
$result/SRR5859757_2.sai \
$dataT/SRR5859757_1.fastq.gz \
$dataT/SRR5859757_2.fastq.gz \
> $result/SRR5859757.sam

time bwa sampe $reference/canFam3.fa \
$result/SRR5867229_1.sai \
$result/SRR5867229_2.sai \
$dataT/SRR5867229_1.fastq.gz \
$dataT/SRR5867229_2.fastq.gz \
> $result/SRR5867229.sam

# Normal
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR5859760_1.fastq.gz > $result/SRR5859760_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR5859760_2.fastq.gz > $result/SRR5859760_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR5859760_1.sai \
$result/SRR5859760_2.sai \
$dataN/SRR5859760_1.fastq.gz \
$dataN/SRR5859760_2.fastq.gz \
> $result/SRR5859760.sam

time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR5867236_1.fastq.gz > $result/SRR5867236_1.sai
time bwa aln -t 4 $reference/canFam3.fa $dataN/SRR5867236_2.fastq.gz > $result/SRR5867236_2.sai

time bwa sampe $reference/canFam3.fa \
$result/SRR5867236_1.sai \
$result/SRR5867236_2.sai \
$dataN/SRR5867236_1.fastq.gz \
$dataN/SRR5867236_2.fastq.gz \
> $result/SRR5867236.sam


####### Convert sam to bam file #######
samtools view -bS@ 4 $result/SRR5859757.sam > $result/SRR5859757.bam
samtools view -bS@ 4 $result/SRR5859760.sam > $result/SRR5859760.bam
samtools view -bS@ 4 $result/SRR5867236.sam > $result/SRR5867236.bam
samtools view -bS@ 4 $result/SRR5867229.sam > $result/SRR5867229.bam

rm $result/SRR5859757.sam $result/SRR5859760.sam $result/SRR5867236.sam $result/SRR5867229.sam
####### Delete #######
#rm $dataN/*.fastq.gz $dataT/*.fastq.gz
#rm $result/*.sai



## merge_bam ##
samtools merge ${result}/SAMN07376273.bam ${result}/SRR5859760.bam ${result}/SRR5867236.bam
samtools merge ${result}/SAMN07376267.bam ${result}/SRR5859757.bam ${result}/SRR5867229.bam  


####### Picard #######
# get header information
samtools view -H ${result}/SAMN07376273.bam > $result/header_N
samtools view -H ${result}/SAMN07376267.bam > $result/header_T

# exclude unmapped reads based on FLAG
samtools view $result/SAMN07376273.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SAMN07376273-cleaned.sam
cat $result/header_N $result/SAMN07376273-cleaned.sam > $result/fooN
mv $result/fooN $result/SAMN07376273-cleaned.sam
samtools view $result/SAMN07376267.bam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > $result/SAMN07376267-cleaned.sam
cat $result/header_T $result/SAMN07376267-cleaned.sam > $result/fooT
mv $result/fooT $result/SAMN07376267-cleaned.sam

# picard sort
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SAMN07376273-cleaned.sam O=$result/SAMN07376273_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups I=$result/SAMN07376267-cleaned.sam O=$result/SAMN07376267_rg_added_sorted.bam SO=coordinate RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

# picard MarkDuplicates
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SAMN07376273_rg_added_sorted.bam O=$result/SAMN07376273_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SAMN07376273-output.metrics REMOVE_SEQUENCING_DUPLICATES=true
time java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar  MarkDuplicates I=$result/SAMN07376267_rg_added_sorted.bam O=$result/SAMN07376267_rg_added_sorted_dedupped_removed.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$result/SAMN07376267-output.metrics REMOVE_SEQUENCING_DUPLICATES=true





####### Remove unneeded files #######
rm $result/SAMN07376273.sam $result/SAMN07376267.sam 
rm $result/header_N $result/header_T
rm $result/SAMN07376273_rg_added_sorted.bam $result/SAMN07376267_rg_added_sorted.bam
rm $result/SAMN07376273-cleaned.sam $result/SAMN07376267-cleaned.sam





####### GATK Realign #######
# Generating interval file for sort.bam
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SAMN07376273_rg_added_sorted_dedupped_removed.bam -o $result/SAMN07376273_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference/canFam3.fa -I $result/SAMN07376267_rg_added_sorted_dedupped_removed.bam -o $result/SAMN07376267_rg_added_sorted_dedupped_removed.bam.intervals --allow_potentially_misencoded_quality_scores

# realign
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SAMN07376273_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SAMN07376273_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores
time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reference/canFam3.fa -I $result/SAMN07376267_rg_added_sorted_dedupped_removed.bam -targetIntervals $result/SAMN07376267_rg_added_sorted_dedupped_removed.bam.intervals -o $result/SAMN07376267_rg_added_sorted_dedupped_removed.realigned.bam --allow_potentially_misencoded_quality_scores





######## Mutect #######
# Current Mutect version uses Java 1.7, while GATK requires Java 1.8. 
module load MuTect/1.1.7-Java-1.7.0_80

time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/dbsnp_9615.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--input_file:normal $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $result/SAMN07376267_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/ND10_441_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/ND10_441_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf

####### Mutect2 #######
#time java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T MuTect2 \
#-R $reference/canFam3.fa \
#-I:tumor $result/ERR1681500_rg_added_sorted_dedupped_removed.realigned.bam \
#-I:normal $result/SAMN07376273_rg_added_sorted_dedupped_removed.realigned \
#--dbsnp $reference/dbsnp_9615.vcf \
#-L $reference/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
#-o $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf





####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar/convert2annovar.pl -format vcf4old $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput

# annovar annotate
perl $reference/annovar/annotate_variation.pl --buildver canFam3 $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $reference






####### remove unneeded files #######
rm $result/*_sorted_dedupped_removed.bam
rm $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS $result/ND10_441_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
rm $result/*.bai
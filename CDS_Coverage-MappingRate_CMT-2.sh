#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N CMT-2
#PBS -l nodes=1:ppn=1
#PBS -l walltime=15:00:00
#PBS -l mem=190gb
#PBS -M kh31516@uga.edu
#PBS -m ae

data='/scratch/kh31516/Mammary/results/CMT-2'
source='/scratch/kh31516/Melanoma/Melanoma_source'

module load SAMtools/1.9-foss-2016b
module load BEDTools/2.26.0-foss-2016b
#module load BEDOPS/2.4.30
module load Subread/1.6.2

cd $data

### Normal samples extrect correctly mapped reads and sort ###
samtools view SRR7780922.bam > SRR7780922.sam
samtools view -H SRR7780922.bam > header_SRR7780922
 
cat SRR7780922.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > SRR7780922-cleaned.sam
cat header_SRR7780922 SRR7780922-cleaned.sam > foo
mv foo SRR7780922-cleaned.sam
samtools view -bS SRR7780922-cleaned.sam > SRR7780922-cleaned.bam
samtools index SRR7780922-cleaned.bam
samtools sort -o SRR7780922_cleaned_sorted.bam SRR7780922-cleaned.bam



### Coverage ###
bedtools coverage -d -a $source/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_uniq_sorted.bed \
-b SRR7780922_cleaned_sorted.bam > SRR7780922_rg_added_sorted.cov

featureCounts -a $source/canFam3.gtf -o SRR7780922_cleaned_sorted.bam.featureCounts \
 -p -B -C --ignoreDup -t exon \
 -Q 1 SRR7780922_cleaned_sorted.bam

rm SRR7780922.sam header_SRR7780922 SRR7780922-cleaned.sam


### Tumor samples extrect correctly mapped reads and sort ###
samtools view SRR7780923.bam > SRR7780923.sam
samtools view -H SRR7780923.bam > header_SRR7780923
 
cat SRR7780923.sam | awk 'BEGIN{FS="\t";OFS="\t"}{if ($2%16==3){print $0}}' > SRR7780923-cleaned.sam
cat header_SRR7780923 SRR7780923-cleaned.sam > foo
mv foo SRR7780923-cleaned.sam
samtools view -bS SRR7780923-cleaned.sam > SRR7780923-cleaned.bam
samtools index SRR7780923-cleaned.bam
samtools sort -o SRR7780923_cleaned_sorted.bam SRR7780923-cleaned.bam


### Coverage ###
bedtools coverage -d -a $source/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list_uniq_sorted.bed \
-b SRR7780923_cleaned_sorted.bam > SRR7780923_rg_added_sorted.cov

featureCounts -a $source/canFam3.gtf -o SRR7780923_cleaned_sorted.bam.featureCounts \
 -p -B -C --ignoreDup -t exon \
 -Q 1 SRR7780923_cleaned_sorted.bam

rm SRR7780923.sam header_SRR7780923 SRR7780923-cleaned.sam




#PBS -S /bin/bash
#PBS -N PRJEB32096-ERR3284983
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=45gb
#PBS -q batch
#PBS -M kh31516@uga.edu
#PBS -m ae


result='/scratch/kh31516/Trio-family/Original_Trio/results/PRJEB32096/' #
reference='/scratch/kh31516/Melanoma/Melanoma_source' 
data='/scratch/kh31516/Trio-family/Original_Trio/data/PRJEB32096/ERR3284983'
annovar_index='/scratch/kh31516/Melanoma_source/annovar'

module load SAMtools/1.9-foss-2016b
module load SRA-Toolkit/2.9.1-centos_linux64
module load BWA/0.7.17-foss-2016b
module load picard/2.16.0-Java-1.8.0_144
module load GATK/3.8-1-Java-1.8.0_144
#module load MuTect/1.1.7-Java-1.7.0_80
module load Perl/5.26.1-GCCcore-6.4.0


####### Download #######
mkdir $result
mkdir $data
#mkdir $data
cd $data
fastq-dump --split-files --gzip ERR3284983
#/scratch/kh31516/Pan_cancer/others/PRJEB32096/scripts/download_sra_sample.sh ERR3284983 $sraRunTable

####### BWA Mapping #######
cd $result
# Tumor
time bwa aln $reference/canFam3.fa $data/ERR3284983_1.fastq.gz > $result/ERR3284983_1.sai
time bwa aln $reference/canFam3.fa $data/ERR3284983_2.fastq.gz > $result/ERR3284983_2.sai


time bwa sampe $reference/canFam3.fa \
$result/ERR3284983_1.sai \
$result/ERR3284983_2.sai \
$data/ERR3284983_1.fastq.gz \
$data/ERR3284983_2.fastq.gz \
> $result/ERR3284983.sam


####### Convert sam to bam file #######

cd $result/


module load SAMtools/0.1.19-foss-2016b

python /scratch/kh31516/Trio-family/Original_Trio/Summarize_BWA_sam.py $result/ERR3284983.sam ERR3284983

cat ERR3284983_Mapping_summary.txt >> $result/PRJEB32096_Mapping_summary.txt

#rm $result/ERR3284983.sam

samtools view -bS $result/ERR3284983.sam > $result/ERR3284983.bam
#!/bin/bash
#PBS -q batch
#PBS -N re-run-WGS-i_8A0A
#PBS -l nodes=1:ppn=4
#PBS -l walltime=300:00:00
#PBS -l mem=60gb
#PBS -M kh31516@uga.edu 
#PBS -m ae



module load Trimmomatic/0.36-Java-1.8.0_144

time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 \
$data/SRR7779460_1.fastq.gz $data/SRR7779460_2.fastq.gz \
SRR7779460_lane1_forward_paired.fq.gz SRR7779460_lane1_forward_unpaired.fq.gz \
SRR7779460_lane1_reverse_paired.fq.gz SRR7779460_lane1_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


cd $SLURM_SUBMIT_DIR
data='/home/kh31516/8940/Illumina_PE_reads'
result='/home/kh31516/8940/results'
cd /home/kh31516/8940/
ml Trimmomatic/0.36-Java-1.8.0_144

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $data/Cbai_G_1.fastq.gz $data/Cbai_G_2.fastq.gz \
$result/Cbai-1_paired.fq $result/Cbai-1_unpaired.fq $result/Cbai-2_paired.fq $result/Cbai-2_unpaired.fq \
ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40 

## change to 4:30
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 $data/Cbai_G_1.fastq.gz $data/Cbai_G_2.fastq.gz \
$result/4-30-Cbai-1_paired.fq $result/4-30-Cbai-1_unpaired.fq $result/4-30-Cbai-2_paired.fq $result/4-30Cbai-2_unpaired.fq \
ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:40 



ml SAMtools/1.6-foss-2016b 
ml BWA/0.7.15-foss-2016b
bwa index Reference_Genome.fasta
bwa mem Reference_Genome.fasta File_1.fastq File_2.fastq > out_algn_mem.sam samtools view –b –S out_algn_mem.sam > out_algn_mem.bam
samtools sort -o out_algn_mem.sorted.bam out_algn_mem.bam
samtools index out_algn_mem.s

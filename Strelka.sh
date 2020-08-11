### Strelka2 github
### https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md

# installation
strelka_folder='/scratch/yf94402/Pan_cancer/source/Strelka'
tumor='/scratch/yf94402/Pan_cancer/result/Mammary/SRR7779460/Trimmomatic_Jan9'
normal='/scratch/yf94402/Pan_cancer/result/Mammary/SRR7779461/Trimmomatic_Jan9'
source='/scratch/yf94402/variant_calling/source/DLA_fasta/DLA_db_genbank_fa/tmp'

module load SAMtools/1.10-GCC-8.2.0-2.31.1

cd $strelka_folder


# download strelka binary
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
# decompress
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
# run demo to check successful installation
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash

# index your fasta file
cd $source
samtools faidx /scratch/kh31516/Melanoma_source/canFam3.fa

# execution
# step1 will generate a "runWorkflow.py" script for you to execute
$strelka_folder/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam $tumor/SRR7779460_DLAdb-sorted.bam \
    --tumorBam $normal/SRR7779461_DLAdb-sorted.bam \
    --referenceFasta $source//scratch/kh31516/Melanoma_source/canFam3.fa \
    --runDir demo_somatic

# your result will be generated in demo_somatic/results/variants
demo_somatic/runWorkflow.py -m local -j 20
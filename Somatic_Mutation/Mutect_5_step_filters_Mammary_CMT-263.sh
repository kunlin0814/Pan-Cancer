#PBS -S /bin/bash
#PBS -q batch
#PBS -N glioma
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -l mem=50gb

data='/scratch/kh31516/Original_Mammary/results/CMT-263'
result='/scratch/kh31516/Pan_cancer/glioma/results/WES'
reference='/work/szlab/Lab_shared_PanCancer/source'
script='/work/szlab/Lab_shared_PanCancer/script'



module load SAMtools/1.9-foss-2016b
module load MuTect/1.1.7-Java-1.7.0_80

#mkdir $result

: "
# Mutect 
time java -jar /usr/local/apps/eb/MuTect/1.1.7-Java-1.7.0_80/mutect-1.1.7.jar --analysis_type MuTect \
--reference_sequence $reference/canFam3.fa \
--dbsnp $reference/SRZ189891_722g.990.SNP.INDEL.CanFam3.1.99.CDS_trimmed.vcf \
--dbsnp $reference/dbSNP_2020/DbSNP_canFam3_version151-DogSD_Feb2020_V3.vcf \
--intervals $reference/Canis_familiaris.CanFam3.1.99.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list \
--defaultBaseQualities 30 \
--input_file:normal $data/SRR7780976_rg_added_sorted_dedupped_removed.realigned.bam \
--input_file:tumor $data/SRR7780979_rg_added_sorted_dedupped_removed.realigned.bam \
--out $result/CMT-263_rg_added_sorted_dedupped_removed.bam_call_stats.txt \
--coverage_file $result/CMT-263_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt \
--vcf $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf
"
####### Annovar #######
# Extract PASS records from vcf
#awk '$7 == "PASS" {print $0}' $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf > $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# 5 Steps filtering
#grep -w KEEP $result/CMT-263_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > $result/CMT-263_PASS.stat
#python $script/Filter_MutectStat_5steps.py $result/CMT-263_PASS.stat $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
#perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
#perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 $result/CMT-263_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $reference/annovar_CanFam3.1.99.gtf

# add gene name

while read line1 line2 line3;
do
    cd /scratch/kh31516/Pan_cancer/glioma/results/WES/${line1}
    python $script/Add_GeneName_N_Signature.py /scratch/kh31516/Pan_cancer/glioma/results/WES/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt
done < /scratch/kh31516/Pan_cancer/glioma/source/PRJNA579792_WES_pairs.txt
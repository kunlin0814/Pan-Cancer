#! /bin/bash

ml Anaconda3/2020.02

# MC
## 5 steps filtering
reference='/work/szlab/Lab_shared_PanCancer/source' 
scripts='/home/kh31516/kh31516_Lab_Share_script'
total_file='/scratch/kh31516/Original_Mammary/source/K9_WXS_mammary_SRR_pairs.txt'
mutect_folder='/scratch/kh31516/Original_Mammary/store/PRJNA489159/newMutectSNP'
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/MT_VAF/Mutect1'
annovar_index='/work/szlab/Lab_shared_PanCancer/source/annovar_CanFam3.1.99.gtf'

while read line1 line2 line3
do
    #echo $line1
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file

#OM
total_file='/scratch/kh31516/Original_Melanoma/source/combined_melanoma.txt'
mutect_folder='/scratch/kh31516/Original_Melanoma/store/PRJEB12081/newMutectSNP'
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/OM_VAF/Mutect1'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


# GLM
total_file='/scratch/kh31516/Pan_cancer/glioma/source/Glioma_WES_Pair.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/glioma/results/store/WES/new_Mutect_SNP'
vaf_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WES/Mutect1VAF'


while read line1 line2 line3
do
     mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file



# HSA
total_file='/scratch/kh31516/Pan_cancer/HSA/source/PRJNA552034-HSA-Normal-Tumor-pair'
mutect_folder='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA552034/new_mutect_results'
vaf_out='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA552034/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
done < $total_file



# LYM
total_file='/scratch/kh31516/Pan_cancer/lymphoma/source/Lym_pairs.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/lymphoma/store/Mutect_Mar25/'
vaf_out='/scratch/kh31516/Pan_cancer/lymphoma/store/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


#OSA
total_file='/scratch/kh31516/Pan_cancer/Osteo/source/Bur_osteo_pair'
mutect_folder='/scratch/kh31516/Pan_cancer/Osteo/store/PRJNA391455/Mutect_Mar25_results'
vaf_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJNA391455/Mutect1VAF'



while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1

done < $total_file


#UCL

total_file='/scratch/kh31516/Pan_cancer/others/source/Unclassified_pairs'
mutect_folder='/scratch/kh31516/Pan_cancer/others/store/Mutect_Mar25'
vaf_out='/scratch/kh31516/Pan_cancer/others/store/Mutect1VAF'



while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file


############################# validation dataset # MC*2, OM*1, OSA*2, HSA *1  ###################
## MC lab

total_file='/scratch/kh31516/Original_Mammary/source/lab_WES_pair.txt'
mutect_folder='/scratch/kh31516/Original_Mammary/store/Lab_WES_MC/Mutect'
vaf_out='/scratch/kh31516/Original_Mammary/store/Lab_WES_MC/Mutect1VAF'

while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file

## MC SU

total_file='/scratch/kh31516/Original_Mammary/source/PRJNA552905_MC_pairs.txt'
mutect_folder='/scratch/kh31516/Original_Mammary/store/PRJNA552905/Mutect'
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA552905/Mutect1VAF'

while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file

## OM sanger

total_file='/scratch/kh31516/Original_Melanoma/source/PRJEB7540_OM.txt'
mutect_folder='/scratch/kh31516/Original_Melanoma/store/PRJEB7540/Mutect'
vaf_out='/scratch/kh31516/Original_Melanoma/store/PRJEB7540/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1

    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file

## OSA
# Tgen

total_file='/scratch/kh31516/Pan_cancer/Osteo/source/Tgen_TrueN_Tpairs_WES.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/Osteo/Tgen_OSA/results_SG'
vaf_out='/scratch/kh31516/Pan_cancer/Osteo/Tgen_OSA/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    echo $line1

    cd $mutect_folder/$line1
    #mv ${line1}*_rg_added_sorted_dedupped_removed.MuTect.vcf ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf
    #mv ${line1}*_PASS.stat ${line1}_PASS.stat
    #mv ${line1}*_rg_added_sorted_dedupped_removed.bam_call_stats.txt ${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt


    #awk '$7 == "PASS" {print $0}' ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
    #perl $annovar_index/convert2annovar.pl -format vcf4old ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
    #perl $annovar_index/annotate_variation.pl --buildver canFam3 ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index
    #python2 $scripts/Add_GeneName_N_Signature.py ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt
    
    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 

    #awk '$7 == "PASS" {print $0}' ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-PASS

    #perl $annovar_index/convert2annovar.pl -format vcf4old ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-PASS > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-PASS-avinput
    #perl $annovar_index/annotate_variation.pl --buildver canFam3 ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-PASS-avinput $annovar_index
    #python2 $scripts/Add_GeneName_N_Signature.py ${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt
    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file
# sanger

total_file='/scratch/kh31516/Pan_cancer/Osteo/source/PRJEB7540_OSA.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Mutect'
vaf_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1
    #awk '$7 == "PASS" {print $0}' ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
    #perl $annovar_index/convert2annovar.pl -format vcf4old ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
    #perl $annovar_index/annotate_variation.pl --buildver canFam3 ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index
    #python2 $scripts/Add_GeneName_N_Signature.py ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt
    
    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file

# HSA 

total_file='/scratch/kh31516/Pan_cancer/HSA/source/PRJNA417727_HSA_WES.txt'
mutect_folder='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA417727/Mutect'
vaf_out='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA417727/Mutect1VAF'


while read line1 line2 line3
do
    mkdir -p $vaf_out/$line1
    cd $mutect_folder/$line1
    #awk '$7 == "PASS" {print $0}' ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS
    #perl $annovar_index/convert2annovar.pl -format vcf4old ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS > ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput
    #perl $annovar_index/annotate_variation.pl --buildver canFam3 ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput $annovar_index
    #python2 $scripts/Add_GeneName_N_Signature.py ${mutect_folder}/${line1}/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function ${reference}/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt
    
    grep -w KEEP $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,4,5,26,27,38,39 > $mutect_folder/$line1/${line1}_PASS.stat

    python $scripts/Filter_MutectStat_5steps.py \
    $mutect_folder/$line1/${line1}_PASS.stat \
    $mutect_folder/$line1/${line1}_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS \
    $mutect_folder/$line1/${line1}_vaf_before.txt \
    $mutect_folder/$line1/${line1}_vaf_after.txt \
    $mutect_folder/$line1/${line1}_whyout.txt \
    $line1 


    mv ${line1}_vaf_before.txt $vaf_out/$line1
    mv ${line1}_vaf_after.txt $vaf_out/$line1
    
done < $total_file


## cat all the file together
final_total_Before_Mutect_output='/scratch/kh31516/Total_VAF_before_Mutect.txt'
final_total_after_Mutect_output='/scratch/kh31516/Total_VAF_after_Mutect.txt'
# MT Kor
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/MT_VAF/Mutect1'
cd $vaf_out

cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="MT" ; print}' >> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="MT" ; print}'>> $final_total_after_Mutect_output

# MT lab


# MT SUN

# OM
# sanger cross
vaf_out='/scratch/kh31516/Original_Mammary/store/PRJNA489159/OM_VAF/Mutect1'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="OM" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="OM" ; print}'>> $final_total_after_Mutect_output

# sanger


# GLM
vaf_out='/scratch/kh31516/Pan_cancer/glioma/results/store/WES/Mutect1VAF'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="GLM" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="GLM" ; print}'>> $final_total_after_Mutect_output
# LYM
vaf_out='/scratch/kh31516/Pan_cancer/lymphoma/store/Mutect1VAF'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="LYM" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="LYM" ; print}'>> $final_total_after_Mutect_output
# HSA
# broad
vaf_out='/scratch/kh31516/Pan_cancer/HSA/store/PRJNA552034/Mutect1VAF'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="HSA" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="HSA" ; print}'>> $final_total_after_Mutect_output

# Upenn

# OSA
# broad
vaf_out='/scratch/kh31516/Pan_cancer/Osteo/store/PRJNA391455/Mutect1VAF'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="OSA" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="OSA" ; print}'>> $final_total_after_Mutect_output

# Tgen


# sanger
/scratch/kh31516/Pan_cancer/Osteo/store/PRJEB7540/Mutect


# UCL
vaf_out='/scratch/kh31516/Pan_cancer/others/store/Mutect1VAF'
cd $vaf_out
cat */*_vaf_before.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="UCL" ; print}'>> $final_total_Before_Mutect_output
cat */*_vaf_after.txt | awk -F"\t" 'BEGIN { OFS = "\t" } {$9="UCL" ; print}'>> $final_total_after_Mutect_output



gzip $final_total_Before_Mutect_output
gzip $final_total_after_Mutect_output
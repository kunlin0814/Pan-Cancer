while read line; 
do sed 's/SRR7779460/'$line'/g' Download_STAR_HISAT2_SRR7779460.sh > scripts/Download_STAR_HISAT2_$line.sh; 
done < /scratch/kh31516/Melanoma/RNA-seq/source/singel_end_total_RNA_seq_melanoma.txt



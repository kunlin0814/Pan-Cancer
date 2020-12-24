#! /bin/bash
while read line1 line2 line3; 
do 

sed 's/i_03A6/'$line1'/g;s/SRR10351810/'$line2'/g;s/SRR10351811/'$line3'/g' 20-WES_whole_Trimmomatic_i_03A6.sh > scripts/20-WES_whole_Trimmomatic_${line1}.sh

done < /scratch/kh31516/Pan_cancer/glioma/source/PRJNA579792_WES_Tumor_NormalPairs+Tumor+normal

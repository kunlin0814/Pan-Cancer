# /bin/bash

while read line ;
do
check=$(grep -w $line /scratch/kh31516/Mammary/results/CMT-*/CMT-*.exonic_variant_function_WithGeneName |cut -f4 |cut -d":" -f5 |sort |uniq -c |sort -nr)
if [ -z "$check" ]; then
    printf "%s\t" "0" >> /scratch/kh31516/Mammary/Mammary_uniq-mut.txt
else
    printf "%s\t" "$check" >> /scratch/kh31516/Mammary/Mammary_uniq-mut.txt
fi
done < /scratch/kh31516/Mammary/source/DNA-repair_list.txt
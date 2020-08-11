#!/bin/bash

############################ description ################################
# This script file accepts a biosample id (SAMNxxxxxxxx), an SRA 
# RunInfo Table, and an (optional) output directory to download all fastq
# files for the given sample. If the output directory isn't supplied, it
# will download the files in the current directory.
# The script will then rename the fastq files to have the biosample id
# instead of the run accession number if the sample has only one run.
# If the sample has more than one run, it will merge the downloaded fastq 
# files and end up with two fastq files named as SAMNxxxxxxxx_1.fastq.gz 
# and SAMNxxxxxxxx_2.fastq.gz and delete the original files.
# Command line:
#       download_sra_sample.sh BiosampleId RunInfoFile [OutputDir]
# Example command line:
#       download_sra_sample.sh SAMN05864126 SraRunTable.txt /scratch/myid/fastq_download
#########################################################################

module load SRA-Toolkit/2.9.1-centos_linux64

echo $(date +"%m-%d-%y %r") Started downloading sample $1

# save the current working directory before changing it
cwd=$(pwd)
if [[ "$#" -eq 3 ]]
then
	# TODO: remove debugging
	echo Changing the current working directory 
	echo "	Running cd $3"
	cd $3
elif [[ "$#" -ne 2 ]]
then
	echo ERROR: Invalid number of arguments passed to download_sra_sample.sh
	exit 1
fi

# no errors and we are now in the right output directory
# get the column number for BioSample and Run
SAMPLE_COL=$(head -1 $2 | tr "\t" "\n" | grep -n "BioSample"'$' | cut -f 1 -d :)
RUN_COL=$(head -1 $2 | tr "\t" "\n" | grep -n "Run"'$' | cut -f 1 -d :)
# get all runs for the sample id
RUN_IDS=($(awk -v sample_id="$1" -v s_col="$SAMPLE_COL" -v r_col="$RUN_COL" '{
  if($s_col == sample_id) {
    print $r_col;
  }
}' FS="\t" "$2"))

# download all runs for this sample using fastq-dump command
echo Downloading sample runs
for i in "${RUN_IDS[@]}"
do
	# TODO: remove debugging
	echo "	Running fastq-dump --split-files --gzip $i"
	fastq-dump --split-files --gzip $i
	echo "	Running rm /scratch/$USER/ncbi/public/sra/${i}.*"
	rm /scratch/$USER/ncbi/public/sra/${i}.*
done

# The final file names
FILE_1_NAME="${1}_1.fastq.gz"
FILE_2_NAME="${1}_2.fastq.gz"

NUMBER_OF_RUNS=len=${#RUN_IDS[@]}
if [[ $NUMBER_OF_RUNS -eq 1 ]]
then
	echo Renaming fastq files from run numbers to sample ids
	# TODO: remove debugging
	echo "	Running mv ${RUN_IDS[0]}_1.fastq.gz $FILE_1_NAME"
	mv ${RUN_IDS[0]}_1.fastq.gz $FILE_1_NAME
	# TODO: remove debugging
	echo "	Running mv ${RUN_IDS[0]}_2.fastq.gz $FILE_2_NAME"
	mv ${RUN_IDS[0]}_2.fastq.gz $FILE_2_NAME
elif [[ $NUMBER_OF_RUNS -gt 1 ]]
then
	# construct input files
	RUN_FASTQ_1=""
	RUN_FASTQ_2=""
	for i in "${RUN_IDS[@]}"
	do
		RUN_FASTQ_1+=" ${i}_1.fastq.gz"
		RUN_FASTQ_2+=" ${i}_2.fastq.gz"
	done
	# run cat, gzip and delete commands
	echo $(date +"%m-%d-%y %r") Merging the run fastq files into ones with sample id and deleting the run files
	# TODO: remove debugging
	echo "	Running zcat$RUN_FASTQ_1 | gzip -c > $FILE_1_NAME"
	zcat$RUN_FASTQ_1 | gzip -c > $FILE_1_NAME
	# TODO: remove debugging
	echo "	Running rm$RUN_FASTQ_1"
	rm$RUN_FASTQ_1
	# TODO: remove debugging
	echo "	Running zcat$RUN_FASTQ_2 | gzip -c > $FILE_2_NAME"
	zcat$RUN_FASTQ_2 | gzip -c > $FILE_2_NAME
	# TODO: remove debugging
	echo "	Running rm$RUN_FASTQ_2"
	rm$RUN_FASTQ_2
else
	echo "ERROR: Couldn't find any runs for sample $1"
	exit 1
fi

# TODO: remove debugging
echo Returning to the current working directory
echo "	Running cd $cwd"
cd $cwd

# The download script has finished
echo $(date +"%m-%d-%y %r") Finished downloading sample $1

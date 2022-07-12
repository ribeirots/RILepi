#!/bin/bash
# Script to split fastq files in blocks and rename them to be uploaded on CHTC to be mapped
# It will split ALL fastq files in the folder. The fastq files CANNOT be compressed.

# It requires Pool Lab script split.pl

# Check if perl script exists.
if [ ! -f split.pl ]
then
	echo "split.pl script not found in this folder."
	echo "Terminating script."
	exit 1
fi

# the following command splits the fastq files in Blocks of 800000 reads per fastq file
perl split.pl 800000

# This step renames the block files to make them compatible with the mapping steps on CHTC.
# This step is very hardwired for the name schemes used in this project, take a look at input file formatting for the CHTC steps.
rename 's/.fast/_R1.fast/' *_R1_*.gz
rename 's/.fast/_R2.fast/' *_R2_*.gz
rename 's/_R1_/_/' *_R1_*.gz
rename 's/_R2_/_/' *_R2_*.gz

ls *R1.fastq.gz > Seq_list_single.txt

sed -i 's/_R1.fastq.gz//' Seq_list_single.txt

# next, upload the Seq_list_single.txt and the block gz files to CHTC. 
# The Seq_list_single.txt file should be on the directory for the mapping.
# The block files should be in a directory could "Seq" within the mapping directory.

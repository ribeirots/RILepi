#!/bin/bash
# pipeline to obtain panel files for both RIL sets
# tribeiro@wisc.edu

# It requires find_snps.pl script
# It requires the _diploid.fasta files from the bam2fasta.sh pipeline
# It requires that the directories RIL_FR and RIL_EF already exist and have the bam file from the RILs

# Check if custom made script exists.
if [ ! -f find_snps.pl ]
then
	echo "find_snps.pl script not found in this folder."
	echo "Terminating script."
	exit 1
fi

# check for the existence of required directories
if [ ! -d RIL_FR/ ]
then
	echo 'RIL_FR does not exist. The RIL bam files should be there as well for the next step.'
	exit
fi

if [ ! -d RIL_EF/ ]
then
	echo 'RIL_EF does not exist. The RIL bam files should be there as well for the next step.'
	exit
fi

## Processing fasta files to generate panel files
for c in X 2L 2R 3L 3R
do
	for l in ZI418N ZI251N EF43N FR320N
	do
		mv ${l}_Chr${c}_diploid.fasta ${l}_Chr${c}.fas
	done
done

# Generate panel files.
for c in X 2L 2R 3L 3R
do
	perl find_snps.pl FR320N_Chr${c}.fas ZI418N_Chr${c}.fas ${c} > RIL_FR/${c}.panel
	perl find_snps.pl EF43N_Chr${c}.fas ZI251N_Chr${c}.fas ${c} > RIL_EF/${c}.panel
done

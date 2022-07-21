#!/bin/bash
# script to generate posterior probability of ancestry with AncestryHMM
# It takes one argument: the RIL set being used: ER or FR

# It requires find_snps.pl script
# It requires rec_map.py script
# It requires recombination map files for all 5 chrm arms (e.g. X_recomb_ER/FR.csv) - only tests for the existence of one of them. Note that Ethiopia and France RILs use different recombination maps
# It requires ../Ancestry_HMM/ 

if [ "$#" -ne 1 ]
then
  echo "The script requires exactly one argument, either ER or FR, representing the RIL set being analyzed."
  exit 1
fi

# Check if required script exists.
if [ ! -f populate_snp_matrix.pl ]
then
	echo "populate_snp_matrix.pl script not found in this folder."
	echo "Terminating script."
	exit 1
fi

if [ ! -f rec_map.py ]
then
	echo "rec_map.py script not found in this folder."
	echo "Terminating script."
	exit 1
fi

if [ ! -f 2L_recomb_ER.csv ] && [ ! -f 2L_recomb_FR.csv ]
then
	echo "Neither 2L_recomb_ER.csv or 2L_recomb_FR.csv files were found in this folder."
	echo "Terminating script."
	exit 1
fi

# check if Ancestry_HMM is present
if [ ! -d ../Ancestry_HMM/ ]
then
	echo "Ancestry_HMM was not found in the directory before (../) this one."
	echo "Terminating script"
	exit 1
fi

# A series of folders should exist to receive the posterior outputs, this block will check whether they exist and will create them if needed.
if [ ! -d posteriors/ ]
then
	mkdir posteriors
fi

if [ ! -d posteriors/X/ ]
then
	mkdir posteriors/X
fi

if [ ! -d posteriors/2L/ ]
then
	mkdir posteriors/2L
fi

if [ ! -d posteriors/2R/ ]
then
	mkdir posteriors/2R
fi

if [ ! -d posteriors/3L/ ]
then
	mkdir posteriors/3L
fi

if [ ! -d posteriors/3R/ ]
then
	mkdir posteriors/3R
fi

# list the RIL bam files
ls *.bam > RIL_list.txt # list all the bam files (all the bam files should be from the RILs to be analized)
sed -i 's/_realign.bam//' RIL_list.txt # edit the file to use it as input to the for loop below

# commands below modified from RIL_pipeline_example.sh, matt's script
for RIL in $(cat RIL_list.txt) #Text file list for each of the bam file prefix names, specify directory and location below
do	
	### Index the line (produces .bam.bai file in the directory containing bams) This is used by mpileup call
	samtools index ${RIL}_realign.bam # not rmdup as sufix anymore, but these files have had duplicates removed already.
	
	#### Generate mpileups. Each number in names corresponds to chromosome region, respectively: X, 2L, 2R, 3L, 3R
	#### See Panel_file_information.docx on FR_ZI_RIL_ancestry shared folder for information on how to generate panel files
	### Generate mpileups
	names="4 3 7 5 8"
	for name in $names
	do
		samtools mpileup -q 20 -r $name ${RIL}_realign.bam > ${name}.mpileup 
	done

	### X Chromosome
	perl populate_snp_matrix.pl X.panel < 4.mpileup > X.ahmm_in.panel
	python3 rec_map.py X_recomb_${1}.csv X.ahmm_in.panel # adjust recombination maps for RILs
	mv X.ahmm_in.panel2 X.ahmm_in.panel

	### 2L
	perl populate_snp_matrix.pl 2L.panel < 3.mpileup > 2L.ahmm_in.panel
	python3 rec_map.py 2L_recomb_${1}.csv 2L.ahmm_in.panel # adjust recombination maps for RILs
	mv 2L.ahmm_in.panel2 2L.ahmm_in.panel

	### 2R
	perl populate_snp_matrix.pl 2R.panel < 7.mpileup > 2R.ahmm_in.panel
	python3 rec_map.py 2R_recomb_${1}.csv 2R.ahmm_in.panel # adjust recombination maps for RILs
	mv 2R.ahmm_in.panel2 2R.ahmm_in.panel

	### 3L
	perl populate_snp_matrix.pl 3L.panel < 5.mpileup > 3L.ahmm_in.panel
	python3 rec_map.py 3L_recomb_${1}.csv 3L.ahmm_in.panel # adjust recombination maps for RILs
	mv 3L.ahmm_in.panel2 3L.ahmm_in.panel

	### 3R
	perl populate_snp_matrix.pl 3R.panel < 8.mpileup > 3R.ahmm_in.panel
	python3 rec_map.py 3R_recomb_${1}.csv 3R.ahmm_in.panel # adjust recombination maps for RILs
	mv 3R.ahmm_in.panel2 3R.ahmm_in.panel

	#### Create the ahmm sample file
	ls ${RIL}_realign.bam | perl -pi -e 's/\n/\t2\n/' > ahmm_in.samples #This command may sometimes produces a warning about -i flag, ignore

	#### Run Ancestry_hmm to generate
	chroms="X 2L 2R 3L 3R"
	for name in $chroms
	do
		#### email Matt if there are any questions about ancestry_hmm parameters (ploidy 2 for all here)
		../Ancestry_HMM/src/ancestry_hmm -i ${name}.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5
		### Move files into a folder "posteriors" and sub-directories (ex. "2L"). The posteriors must have the naming scheme "<line>.posterior
		### and be in a subdirectory with names X, 2L, 2R, 3L, and 3R for the script used in the next steps.
		mv ${RIL}_realign.bam.posterior posteriors/${name}/${RIL}.posterior 
	done
	
	### Remove unecessary intermediate files
	rm 3.mpileup
	rm 4.mpileup
	rm 5.mpileup
	rm 7.mpileup
	rm 8.mpileup
	rm *.ahmm_in.panel
	rm ahmm_in.samples
done

#### Following this, place file RIL_genotypes_windows.pl in the posteriors folder. Ensure the script is updated at the "@ my RILS"
#### line to include the lines you are using (and are the names of the posterior files. "my @chrs" should be: "X", "2L", "2R", "3L" and "3R",
#### the name of the sub-directories containing ahmm posteriors. Place window files into the posterior folder and run the perl script
#### to generate ancestry calls for each line in a single csv file.

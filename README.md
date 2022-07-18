# Invetigating Adaptive Epistasis using *Dmel* RILs

This is a collaborative project involving genotyping and phenotyping of two sets of *Drosophila melanogaster* RILs. One set consists of a France *vs* Zambia (France RILs) cross and the other an Ethiopia *vs* Zambia (Ethiopia RILs) cross. For each set, two phenotypes have been measured. For the France RILs, the phenotypes were ethanol tolerance and cold tolerance. For the Ethiopia RILs, the phenotypes were wing length and pigmentation. The pigmentation phenotype is further subdivided into three distinct measures: pigmentation score in gray scale for the background of the mesopleura (a thorax trait), and two abdominal traits, a pigmentation score the 4th abdominal segment (A4 Background), and the proportion of the 4th abdominal background covered by a black stripe.

The genomes for each RIL and the parental strains were obtained with Next Gen sequencing using Illumina platform.

## Short read processing
*Scripts for this section are located in the Alignment directory.*

This section contains the steps taken to analyse the genomic data, starting with the fastq files output from the Illumina sequecing.

### DNA Alignment
We used the UW-Madison Center for High Throughput Computing (CHTC) to perform the DNA alignment of our samples. Using CHTC allowed us to greatly parallelize the alignment step. It required that our initial fastq files be split in shorter files and each of these shorter fastq files were submitted as independent jobs to be aligned - and then latter merged back together. The *fastq_processing1.sh* handles the splitting step for each sample and generates a file listing all the shorter fastq blocks which will be used for the parallel submission. After running this script locally, in the same folder with the fastq files to be split, the block files and the newly generated *Seq_list_single.txt* file need to be uploaded to CHTC.

On CHTC, after the files are placed in their appropriate directories, we submitted the jobs with the *align_pt1.sub* script. This script uses as executable file the *align_pt1.sh* script. The *align_pt1.sh* script will perform the DNA alignment using bwa mem. The complete pipeline used in each parallel job also included steps using *samtools* to process the bam files. It uses *collate* to organize the reads, *fixmate* to correct any flaws in read-pairing, *sort* to sort the reads and finally *markdup* with *-r* flag to identify and remove duplicate reads. In the end, we have bam files sorted and without duplicates for each fastq block. The block files need to be download to continue the pipeline locally.

The next step uses a simple script, *block2rmdupbam.sh* to merge the blocks and remove duplicates again if any was still left after the merging. Starting from here, the pipeline uses one script if you are handling the parental genomes and another if you analyzing the RIL genomes. They both start by using GATK to generate a bam files realigning the reads around indels. But while the parental genome script continues to generate fasta files for each *Dmel* chromosome arm (X, 2L, 2R, 3L, and 3R) the script for the RILs stops at the realigned bam files.

To generate the fasta files from the vcf files I am having to create a new pipeline, diverging again from the NEXUS pipuline. Previous, in the NEXUS, the custom made 'VCF\_to\_Seq' perl script would adjust the size of the contigs based on the indel and sites (shifted) vcf. Because my genome sizes were not shifted, since I did not use two rounds of DNA alignment, using indel and sites vcf was resulting in fasta of the wrong size. Now, alternatively, I am using a new custom made script, *indel\_vcf.py*, to read the indel vcf file from GATK and create a vcf file with the location of all the indels - adding 3 extra sites around them - to be filtered out. Then, using _bedtools intersect_, I generate a sites vcf without the indel regions. This file is used as input for _bcftools filter_ to remove low quality (<32, as per criterium of the VCF\_to\_Seq perl script, also remove GT=./. - missing data) sites. Following that, the filtered vcf file will be used to generate consensus sequence with _bcftool consensus_. The output fasta is then split into one fasta per chromosome arm with the script _split\_fasta.py_. 

The RILs realign bamfiles and the fasta files will be used as input for the steps. First, the fasta files will be used to generate ancestry panel files. Then, the RIL realigned bam files and the panel files will be used to generate ancestry calls for each RIL. Let's discuss the steps to generate panel first.

### Ancestry panel files
test2

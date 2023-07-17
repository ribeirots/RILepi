# Invetigating Adaptive Epistasis using *Dmel* RILs

This is a collaborative project involving genotyping and phenotyping of two sets of *Drosophila melanogaster* RILs. One set consists of a France *vs* Zambia (France RILs) cross and the other an Ethiopia *vs* Zambia (Ethiopia RILs) cross. For each set, two phenotypes have been measured. For the France RILs, the phenotypes were ethanol tolerance and cold tolerance. For the Ethiopia RILs, the phenotypes were wing length and pigmentation. The pigmentation phenotype is further subdivided into three distinct measures: pigmentation score in gray scale for the background of the mesopleura (a thorax trait), and two abdominal traits, a pigmentation score the 4th abdominal segment (A4 Background), and the proportion of the 4th abdominal background covered by a black stripe.

The genomes for each RIL and the parental strains were obtained with Next Gen sequencing using Illumina platform.

## DNA Alignment
*Scripts for this section are located in the _Alignment_ directory.*

This section contains the steps taken to analyse the genomic data, starting with the fastq files output from the Illumina sequecing.

We used the UW-Madison Center for High Throughput Computing (CHTC) to perform the DNA alignment of our samples. Using CHTC allowed us to greatly parallelize the alignment step. It required that our initial fastq files be split in shorter files and each of these shorter fastq files were submitted as independent jobs to be aligned - and then later merged back together. The *fastq_processing1.sh* handles the splitting step for each sample and generates a file listing all the shorter fastq blocks which will be used for the parallel submission. After running this script locally, in the same folder with the fastq files to be split, the block files and the newly generated *Seq_list_single.txt* file need to be uploaded to CHTC.

On CHTC, after the files are placed in their appropriate directories, we submitted the jobs with the *align_pt1.sub* script. This script uses as executable file the *align_pt1.sh* script. The *align_pt1.sh* script will perform the DNA alignment using bwa mem. In the end, we have bam files each fastq block. The block files need to be download to continue the pipeline locally.

Locally, the *block2realign.sh* script will be used to merge the blocks and remove duplicated reads. It uses *samtools* to process the bam files. It uses *collate* to organize the reads, *fixmate* to correct any flaws in read-pairing, *sort* to sort the reads and finally *markdup* with *-r* flag to identify and remove duplicate reads. Then, it uses *Picard* and *GATK* (with specific versions and requirements, see script) to re-align the bam file around indels - it indexes the bam file in the process. The script *realign2vcf.sh* is then used to generate the sites vcf and indels vcf files - this could be integrated in the previous script if many bam files will have vcf called, which is the case in this project but not in the SIBSAM project.

The final bam file can be used on the RILs downstream ancestry analysis. The final vcf file can be used for downstream parental genome analysis. 

## Ancestry calls
*Scripts for this section can be found in the _Ancestry_ directory.*

To obtain panel files, the parental vcfs need to be generated first. They are generated using the script _realign2vcf.sh_. Then, the script _vcf2panel.py_ is used, as in the command examples shown in _panel\_cmds.txt_. The _vcf2panel.py_ requires a set of functions present in three other python scripts: _infordepth\_functions.py_, _parentalVCF\_dict.py_, _readPairs.py_.

Order of the genomes matter and herein EF and FR were used first and will have genotype 0. The panel files are stored in their RIL's respective folders along the bam files for all the RILs (already without duplicates and realigned around indels). Then, we call the script _ancestry.sh_ to use the panels and bam files to call AncestryHMM and obtain posterior ancestry probabilities for each RIL. The _ancestry.sh_ script loops through each RIL and each chromosome arm to obtain ril-arm-specific mpileup files (obtained with _samtools mpileup_) and the ril-arm-specific panel files (obtained with the _populate\_snp\_matrix.pl_ script, script that was obtained with the Ancestry HMM pipeline). To adjust the distance between the SNPs accordingly to the RIL experimental design we used the script _rec\_map.py_ and _recomb.csv_ files generated for each RIL set and chromosome arm based on how many generations of recombination they had and on the recombination rates for _Dmel_ from Comeron et al. 2012 (FR RILs are F12 and EF RILs are F13. Therefore, recombination rates were multiplied by 11 and 12, respectively). With the panel files per RIL per arm we call _ancestry\_hmm_ to obtain the posterior probabilities.

Following this, we used the script _RIL\_genotypes\_windows.pl_, placed in the _posteriors_ folder, to call ancestries for windows and not per SNP. The files delimiting the window sizes need to be located in the _posterior_ folders as well. Here, we used _windows\_ZI1k\_Chr*.txt_, which has windows defined to contain 1k SNPs in the ZI population. We have to manually update the genotypes perl script to ensure the script has the correct RILs at the _"@ my RILS"_ variable, we also need to manually input the output file name.

We obtained the ancestry calls in the files: _*\_RIL\_genotypes\_winZI1k.txt_.

## G x P input files
*Scripts for this section can be found in the _GenoPheno_ directory.*

We combined the ancestry genotypes per window obtained above with phenotypic data to generate files that will be analyzed in our mapping and epistasis analyses. We've switched our focus from r/qtl R package to a custom, modified, python implementation of their models. In our implementation, we consider the phenotype of the intermediary genotypes to be the average of the homozygous genotypes (additive model). We used a set of files and scripts to generate the input files, as described below.

We used _input\_rqtl\_nopheno.py_ to transform the data into the input format for rqtl. Before adding phenotypes, we used a script called _rmLinesRqtl.py_ to remove desired lines based on quality control - this script needs to be modified case by case - then used _input\_rqtl\_addpheno.py_ to add the desired phenotype (see next paragraph). We used the script _missing\_genotypes.py_ to imputate some of the missing windows (if they were surrounded by the same genotype, and there were no more than 10 contiguous missing windows, the genotyped was replaced with the surrounding genotypes). We used the script _window\_merge\_genotype.py_ to merge neighboring windows if the genotype for these windows were the same for all the RILs. Then, we used the script _comeron\_rqtl.py_ to add recombination data to the file, basically adding the cM position for the midpoint of each window. We also added to this github directory an example bash script with a pipeline from the ancestry call output from the previous step all the way to the input file for the next analyses, it is called _genopheno\_pipeline.sh_.

The format of the phenotype data we use includes a header in the first row and the RILs and phenotypes in the remaning rows. There are two columns, separated by comma (,). The first column has the RIL ID number (without any letters) and the second column has the phenotype for that RIL. Example:

```
RIL,pheno
3,0.9
5,0.32
11,0.3697
15,954 
```

## Genome Scan - One Window at a time (scanone)
*Scripts for this section can be found in the _rqtl_ directory.*
The output of the previous section is a file containing phenotype and genotype information. It will be the input file for the _scanone.py_ script, our implementation of the homonymous function from r/qtl.  It requires a set of python packages and can only be parallelized "manually" - submitting many repeated jobs changing the "batch" number to obtain unique output files.

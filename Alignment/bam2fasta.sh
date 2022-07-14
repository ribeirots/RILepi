#!/bin/bash
# Script for the parental genomes.
# Script starting with bam files with duplicated removed. The file has to be labeled FILEID_rmdup.bam
# FILEID is the input argument that the script takes
# The final output contains several vcf files and fasta files per chromosome arm for Dmel (X, 2L, 2R, 3L, and 3R)

# It requires picard and GATK, specific versions can be seen in the commands below.
# It requires Dmel reference fasta
# It requires Pool Lab script VCF_to_Seq_diploid_ambiguities.pl
# It requires indel_bed.py script (made by me)
# It requires filter_indel_vcf.py script (made by me)


# Check if custom made script exists.
if [ ! -f VCF_to_Seq_diploid_ambiguities.pl ]
then
	echo "VCF_to_Seq_diploid_ambiguities.pl script not found in this folder."
	echo "Terminating script."
	exit 1
fi

if [ ! -f indel_bed.py ]
then
	echo "indel_bed.py script not found in this folder."
	echo "Terminating script."
	exit 1
fi

if [ ! -f filter_indel_vcf.py ]
then
	echo "filter_indel_vcf.py script not found in this folder."
	echo "Terminating script."
	exit 1
fi

# creates tmp/ directory if one does not already exists.
if [ ! -d tmp/ ]
then
	mkdir tmp
fi

# creates output/ directory if one does not already exists.
if [ ! -d output/ ]
then
	mkdir output
fi

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/picard-tools-1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$1 INPUT=${1}_rmdup.bam OUTPUT=${1}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${1}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ./alignment_software/dmel_ref/DmelRef.fasta -I ${1}_header.bam -o ${1}.intervals 2> ${1}_realigntarget.log

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R ./alignment_software/dmel_ref/DmelRef.fasta -targetIntervals ${1}.intervals -I ${1}_header.bam -o ${1}_realign.bam 2> ${1}_indelrealign.log

rm ${1}_header.bam
rm ${1}.intervals

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ./alignment_software/dmel_ref/DmelRef.fasta -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -I ${1}_realign.bam -o  ${1}_SNPs.vcf

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ./alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ./alignment_software/dmel_ref/DmelRef.fasta -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -minIndelFrac 0.51 -minIndelCnt 3 -glm INDEL -I ${1}_realign.bam -o ${1}_INDELS.vcf

python3 indel_bed.py ${1}_INDELS.vcf

java -Xmx5g -jar ./alignment_software/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${1}_realign.bam

java -Xmx5g -jar ./alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ./alignment_software/dmel_ref/DmelRef.fasta -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -out_mode EMIT_ALL_SITES -I ${1}_realign.bam -o ${1}_sites.vcf 

python3 filter_indel_vcf.py ${1}_INDELS_indels.bed ${1}_sites.vcf


bgzip ${1}_sites.vcf

tabix -p vcf ${1}_sites.vcf.gz

bcftools consensus -H I -a N -f -f ./alignment_software/dmel_ref/DmelRef.fasta ${1}_sites_indelfree.vcf.gz > ${1}_sites.fasta # this fasta file will not be separated by chrm arm, will need further processing. The VCF was filtered from indels and low quality sites. -H I will output IUPAC code from the genotypes and -a N will output N when the site is absent.

#grep -v '#' ${1}_sites.vcf > ${1}_shifted.vcf # this is only called shifted because that is the name format used on the perl script below (previously used in the NEXUS pipeline)
#gzip ${1}_shifted.vcf
#gzip ${1}_sites.vcf
mv ${1}_site*.vcf.g* output/
#perl VCF_to_Seq_diploid_ambiguities.pl

mv ${1}_realign.ba* output/
mv *.fas output/
mv ${1}_INDEL*.vc* output/
mv ${1}_SNPs.vc* output/
mv *.log output/
#mv ${1}_shifted.vcf.g* output/

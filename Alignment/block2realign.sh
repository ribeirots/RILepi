#!/bin/bash
# 1: bamfile (everything except ".bam")
# command example: bash ./block2realign.sh bampool

echo 'requires java8 (on marula: conda activate tiagojav8)'

id=$1

dmel_ref=/raid10/Tiago/PigmCages/scripts/alignment_software/dmel_ref/DmelRef.fasta # Drosophila melanogaster reference genome
align_pack=/raid10/Tiago/PigmCages/scripts/alignment_software # path to Picard and GATK, note that that GATK version used here required java 8. 


# Merge the block files, assumes that they end in "indexed.bam" as per the previous pipeline step
samtools merge ${id}.bam *indexed.bam

# rmdup
samtools collate -o ${id}_bam_collate.bam ${id}.bam
samtools fixmate -m ${id}_bam_collate.bam ${id}_bam_fixmate.bam
samtools sort -o ${id}_bam_position.bam ${id}_bam_fixmate.bam
samtools markdup -r -d 2500 --threads 20 -f stats.txt ${id}_bam_position.bam ${id}_rmdup.bam

rm ${id}_bam_collate.bam
rm ${id}_bam_position.bam
rm ${id}_bam_fixmate.bam

java -Xmx5g -jar $align_pack/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${id}_rmdup.bam

java -Xmx5g -jar $align_pack/picard-tools-1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$id INPUT=${id}_rmdup.bam OUTPUT=${id}_header.bam

java -Xmx5g -jar $align_pack/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${id}_header.bam

java -Xmx5g -jar $align_pack/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $dmel_ref -I ${id}_header.bam -o ${id}.intervals 2> ${id}_realigntarget.log

java -Xmx5g -jar $align_pack/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R $dmel_ref -targetIntervals ${id}.intervals -I ${id}_header.bam -o ${id}_realign.bam 2> ${id}_indelrealign.log

# For parental genomes: proceed to generate panel files
# For offspring: the bam files are the input for the ancestry hmm pipeline


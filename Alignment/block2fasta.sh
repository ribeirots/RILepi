#!/bin/bash
# Only argument, the FILEID before FILEID_Block*.bam
# need to run on java 8. On marulage: conda activate tiagojav8

# this script hasnt been ran from beginning to end in a single click yet.
id=$1

script_path=/raid10/Tiago/RILS/ParentalGenomes/block2fasta_scripts # path to custom scripts used to process the fasta and vcf files to generate the final diploid.fasta output
dmel_ref=/raid10/Tiago/PigmCages/scripts/alignment_software/dmel_ref/DmelRef.fasta # Drosophila melanogaster reference genome
align_pack=/raid10/Tiago/PigmCages/scripts/alignment_software # path to Picard and GATK, note that that GATK version used here required java 8. 
echo "java 8 is needed for the version of GATK used in this script."

if [ ! -f ${id}_Block00001_*.bam ]
then
		echo "block files not found."
			echo "Terminating script."
				exit 1
fi

mkdir ${id}_intermediary_files
tmp_f=${id}_intermediary_files

samtools merge ${id}_merged.bam ${id}_Block*_align_indexed.bam

ls -l *.bam >> log_block2fasta.txt
echo 'done with merge.' > log_block2fasta.txt


### Remove duplicate sequences
samtools collate -o bam_collate.bam ${id}_merged.bam
samtools fixmate -m bam_collate.bam bam_fixmate.bam
samtools sort -o bam_position.bam bam_fixmate.bam
samtools markdup -r -d 2500 -f stats.txt --threads 20 bam_position.bam ${id}_rmdup.bam

ls -l *.bam >> log_block2fasta.txt
echo 'done with rmdup.' >> log_block2fasta.txt
mv *Block*.bam $tmp_f
mv *_merged.bam $tmp_f

if [ ! -f ${id}_rmdup.bam ]
then
		echo "rmdup file not found."
			echo "Terminating script."
				exit 1
fi

java -Xmx5g -jar ${align_pack}/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${id}_rmdup.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ${align_pack}/picard-tools-1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$id INPUT=${id}_rmdup.bam OUTPUT=${id}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ${align_pack}/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${id}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $dmel_ref -I ${id}_header.bam -o ${id}.intervals 2> ${id}_realigntarget.log

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R $dmel_ref -targetIntervals ${id}.intervals -I ${id}_header.bam -o ${id}_realign.bam 2> ${id}_indelrealign.log

ls -l *.bam >> log_block2fasta.txt
echo 'done with indel realignment.' >> log_block2fasta.txt

mv *_rmdup.bam $tmp_f
mv *_header.bam $tmp_f
mv bam_collate.bam $tmp_f
mv bam_fixmate.bam $tmp_f
mv bam_position.bam $tmp_f
mv *.intervals $tmp_f

# realign to vcf 
java -Xmx5g -jar ${align_pack}/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${id}_realign.bam

java -Xmx5g -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $dmel_ref -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -minIndelFrac 0.51 -minIndelCnt 3 -glm INDEL -I ${id}_realign.bam -o ${id}_INDELS.vcf

java -Xmx5g -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $dmel_ref -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -I ${id}_realign.bam -o ${id}_SNPs.vcf

java -Xmx5g -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $dmel_ref -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -out_mode EMIT_ALL_SITES -I ${id}_realign.bam -o ${id}_sites.vcf

ls -l *.vcf >> log_block2fasta.txt
echo 'Done generating 3 vcf files with GATK.' >> log_block2fasta.txt

python3 ${script_path}/indel_vcf.py ${id}_INDELS.vcf # generates vcf file with all the sites to be filtered out for being indels or near indels (3 bp). The output is used in the bedtools intersect command below

echo 'Done running indel_vcf.py.' >> log_block2fasta.txt

bedtools intersect -v -a ${id}_sites.vcf -b ${id}_INDELS_indelfilter.vcf > ${id}_noindel_sites.vcf

ls -l *_noindel_sites.vcf >> log_block2fasta.txt
echo 'Done generating noindel_sites.vcf.' >> log_block2fasta.txt

bgzip ${id}_noindel_sites.vcf
tabix -p vcf ${id}_noindel_sites.vcf.gz

bgzip ${id}_sites.vcf
tabix -p vcf ${id}_sites.vcf.gz

bcftools filter -i 'QUAL<32 | GT="./."' ${id}_sites.vcf.gz -o ${id}_to_exclude.vcf # marks the regions to be filtered out (no data and low qual)
ls -l *_to_exclude.vcf >> log_block2fasta.txt
echo 'Done generating _to_exclude.vcf.' >> log_block2fasta.txt

bgzip ${id}_to_exclude.vcf
bgzip ${id}_INDELS_indelfilter.vcf # output of indel_vcf.py

tabix -p vcf ${id}_to_exclude.vcf.gz
tabix -p vcf ${id}_INDELS_indelfilter.vcf.gz

bcftools consensus -I -m ${id}_to_exclude.vcf.gz -f /raid10/Tiago/PigmCages/scripts/alignment_software/dmel_ref/DmelRef.fasta ${id}_sites.vcf.gz > ${id}_sites.fasta  # this fasta file will not be separated by chrm arm, will need further processing. The VCF was filtered from indels. This command should filter out low quality sites (-i should only include QUAL>=32). -I will output IUPAC code from the genotypes. The '-a N', that used to output N when the site is absent, is obsolete now. I will need to check the output to see what was done with the absent sites. Manual for the 1.8 version of bcftools: http://www.htslib.org/doc/1.8/bcftools.html#common_options

ls -l *.fasta >> log_block2fasta.txt
echo 'Done generating site.fasta.' >> log_block2fasta.txt

python3 ${script_path}/split_fasta.py ${id}_sites.fasta ${id}

ls -l *.fasta >> log_block2fasta.txt
echo 'Done splitting fasta.' >> log_block2fasta.txt

gzip -d ${id}_sites.vcf.gz

python3 ${script_path}/iupac_fasta.py ${id} X
python3 ${script_path}/iupac_fasta.py ${id} 2L
python3 ${script_path}/iupac_fasta.py ${id} 2R
python3 ${script_path}/iupac_fasta.py ${id} 3L
python3 ${script_path}/iupac_fasta.py ${id} 3R

ls *.* >> log_block2fasta.txt

mv *_INDELS.vcf $tmp_f
mv *_SNPs.vcf $tmp_f
mv *.vcf.gz $tmp_f

bgzip ${id}_sites.vcf

mv *ChrX.fasta $tmp_f
mv *Chr2L.fasta $tmp_f
mv *Chr2R.fasta $tmp_f
mv *Chr3L.fasta $tmp_f
mv *Chr3R.fasta $tmp_f
mv *sites.fasta $tmp_f

mv ${id}_ChrX_diploid.fasta ${id}_ChrX_diploid.fas
mv ${id}_Chr2L_diploid.fasta ${id}_Chr2L_diploid.fas
mv ${id}_Chr2R_diploid.fasta ${id}_Chr2R_diploid.fas
mv ${id}_Chr3L_diploid.fasta ${id}_Chr3L_diploid.fas
mv ${id}_Chr3R_diploid.fasta ${id}_Chr3R_diploid.fas

mv *.idx $tmp_f
mv *.tbi $tmp_f
mv *.bai $tmp_f
mv *.log $tmp_f

rm -r tmp

echo 'Script finished.' >> log_block2fasta.txt


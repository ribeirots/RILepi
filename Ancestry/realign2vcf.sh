#!/bin/bash
align_pack=/raid10/Tiago/PigmCages/scripts/alignment_software # path to Picard and GATK, note that that GATK version used here required java 8. 
dmel_ref=/raid10/Tiago/PigmCages/scripts/alignment_software/dmel_ref/DmelRef.fasta # Drosophila melanogaster reference genome

id=$1

java -Xmx5g -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $dmel_ref -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -out_mode EMIT_ALL_SITES -I ${id}_realign.bam -o ${id}_sites.vcf


java -Xmx5g -jar ${align_pack}/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $dmel_ref -mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 2 -minIndelFrac 0.51 -minIndelCnt 3 -glm INDEL -I ${id}_realign.bam -o ${id}_INDELS.vcf

#!/bin/bash

for SAMPLE in TWM_17_346
do

echo 'Performing Mutect2 on' ${SAMPLE}

gatk Mutect2 -I RNA_Exomes_1_27_20/${SAMPLE}_bqsr.bam -R mutation_working_files/human_g1k_v37.fasta -O ${SAMPLE}_mutect.vcf --f1r2-tar-gz ${SAMPLE}_f1r2.tar.gz --max-mnp-distance 0 

gatk LearnReadOrientationModel -I ${SAMPLE}_f1r2.tar.gz -O ${SAMPLE}_read-orientation-model.tar.gz

gatk GetPileupSummaries -I ${SAMPLE}_bqsr.bam -V mutation_working_files/somatic-b37_small_exac_common_3.vcf -L mutation_working_files/somatic-b37_small_exac_common_3.vcf -O ${SAMPLE}_getpileupsummaries.table

gatk CalculateContamination -I ${SAMPLE}_getpileupsummaries.table -O ${SAMPLE}_calculatecontamination.table

gatk FilterMutectCalls -V ${SAMPLE}_mutect.vcf --contamination-table ${SAMPLE}_calculatecontamination.table --ob-priors ${SAMPLE}_read-orientation-model.tar.gz -O ${SAMPLE}_mutect_filtered.vcf -R mutation_working_files/human_g1k_v37.fasta --max-alt-allele-count 1

snpeff -stats ${SAMPLE}.html -v GRCh37.75 ${SAMPLE}_mutect_filtered.vcf > ${SAMPLE}_mutect_filtered_b37_ann.vcf

vcf2tsv -g ${SAMPLE}_mutect_filtered_b37_ann.vcf > ${SAMPLE}_mutect_b37_ann.tsv

done
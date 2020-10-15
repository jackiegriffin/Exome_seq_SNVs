#!/bin/bash

#!/bin/bash

for SAMPLE in TWM_17_049 TWM_17_346 TWM_17_347 TWM_17_348

do
    echo 'Performing Base Recalibration on' ${SAMPLE}

    gatk BaseRecalibrator -I ${SAMPLE}_rmdup.bam -R ../mutation_working_files/human_g1k_v37.fasta -O ${SAMPLE}_recal_data.table --known-sites ../mutation_working_files/dbsnp_138.b37.vcf

    gatk ApplyBQSR -I ${SAMPLE}_rmdup.bam -R ../mutation_working_files/human_g1k_v37.fasta --bqsr-recal-file ${SAMPLE}_recal_data.table -O ${SAMPLE}_bqsr.bam

done
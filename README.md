# Exome-mutation-calls

Analyze exome seq data to identify SNVs, InDels, and CNAs.

Interested in determining:
1. Genetic variability among the 3 tumors at baseline
2. Genetic changes observed on Day 90 vs. baseline
3. Genetic changes observed in estrogen-independent tumors vs. baseline
 
There are 11 tumor samples, all based on MCF-7 cells, grouped into 3 time points:

TWM-17-049: MCF-7 tumor growing on E2 (baseline)
TWM-17-050: MCF-7 tumor growing on E2 (baseline)
TWM-17-051: MCF-7 tumor growing on E2 (baseline)
TWM-17-238: MCF-7 tumor following 90 d of E2 withdrawal
TWM-17-240: MCF-7 tumor following 90 d of E2 withdrawal
TWM-17-241: MCF-7 tumor following 90 d of E2 withdrawal
TWM-18-293: MCF-7 tumor that regrew without E2 (estrogen-independent)
TWM-18-294: MCF-7 tumor that regrew without E2 (estrogen-independent)
TWM-18-346: MCF-7 tumor that regrew without E2 (estrogen-independent)
TWM-18-347: MCF-7 tumor that regrew without E2 (estrogen-independent)
TWM-18-348: MCF-7 tumor that regrew without E2 (estrogen-independent)


Call variants and indicate:
- gene
- allelic frequency
- chromosome location
- nucleotide change
- amino acid change

***********************************************************************************
 
 Previously analyzed by Jason, raw fastq files were:
 
 1) sorted
 2) PCR duplicates removed
 3) BWA-MEM aligned
 
 ***********************************************************************************

I received BWA-MEM aligned bqsr.bam files to call mutations in remaining half os samples:

TWM-17-049
TWM-17-050
TWM-17-051
TWM-17-238
TWM-17-240
TWM-17-241
TWM-17-293
TWM-17-294
TWM-17-346
TWM-17-347
TWM-17-348


STEPS:

Install software needed to run bashscript:

1. $ brew install gatk
2. vcf2tsv through github (https://github.com/sigven/vcf2tsv)
   vcf2tsv dependencies:
    - Python 3.7.6
        $ python3 --version
    - cyvcf2
        $ pip install cyvcf2
    - numpy (https://scipy.org/install.html)
        $ conda install -c anaconda numpy
 3. snpeff through sourceforge (http://snpeff.sourceforge.net/download.html)


***************************************************************************************

gatk Mutec2 bash script ---- w/o vcf2 and snpeff lines ----

***************************************************************************************

#!/bin/bash

for SAMPLE in TWM_17_049
do

echo 'Performing Mutect2 on' ${SAMPLE}

gatk Mutect2 -I RNA_Exomes_1_27_20/${SAMPLE}_bqsr.bam -R mutation_working_files/human_g1k_v37.fasta -O ${SAMPLE}_mutect.vcf --f1r2-tar-gz ${SAMPLE}_f1r2.tar.gz --max-mnp-distance 0 

gatk LearnReadOrientationModel -I ${SAMPLE}_f1r2.tar.gz -O ${SAMPLE}_read-orientation-model.tar.gz

gatk GetPileupSummaries -I ${SAMPLE}_bqsr.bam -V mutation_working_files/somatic-b37_small_exac_common_3.vcf -L mutation_working_files/somatic-b37_small_exac_common_3.vcf -O ${SAMPLE}_getpileupsummaries.table

gatk CalculateContamination -I ${SAMPLE}_getpileupsummaries.table -O ${SAMPLE}_calculatecontamination.table

gatk FilterMutectCalls -V ${SAMPLE}_mutect.vcf --contamination-table ${SAMPLE}_calculatecontamination.table --ob-priors ${SAMPLE}_read-orientation-model.tar.gz -O ${SAMPLE}_mutect_filtered.vcf -R mutation_working_files/human_g1k_v37.fasta --max-alt-allele-count 1

done

***************************************************************************************














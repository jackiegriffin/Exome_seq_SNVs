# Exome-mutation-calls

Please analyze exome seq data to identify SNVs, InDels, and CNAs.

Todd is most interested in determining:
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
    - python 3
    - cyvcf2
    - numpy
    
    
 
sudo apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev curl



















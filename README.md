# Exome-mutation-calls

Analyze exome seq data to identify SNVs, InDels, and CNAs.

Interested in..
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

I received BWA-MEM aligned bqsr.bam files to call mutations for samples:

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


***********************************************************************************

RUN gatk-Mutec2 bash script 'RNA_Exomes.sh'

************************************************************************************

PROCESS VARIANTs - call using 'mutect_process' function in R

    - 'hgvs_p' column = protein changes. REMOVE ROWS WITH NO PROTEIN CHANGES; cells containing $""$
    - 'filter' column = only keep cells labeled $'PASS'$ [filter = pass keep]
    - remove synonymous variants
    - AF = mutant allelic frequency
    - DR = read depth

***********************************************************************************

CONSIDERATIONS:

- What is the biologic relevance of called variants ?
- Look at mutation burden and how they compare to data in similar diseases
- Show nucleotide changes are relevant to other similar tumor related data.

- Variant allele frequency greater than 40% might be a germline mutation. so lesss thatn 45 % is condisered a somatic mutation. threshold of elow 1% could be noise. 

Bonnie laus threshold to call variant calls was between 0.1% and 10

***

identify driver mutations through germlne mutation (greater than 40$ allelic frequency of mutational brden from exome seq
vs
identify passenger genes/genomic instability through somatic mutation (between 0.1-10%)


aytologous samples from patients for engraftment into mice for murine models






















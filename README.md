# Exome-mutation-calls


# Upload data 2----
filtered_variant_calls <- list.files(path = "C:/Users/jackie/OneDrive - Dartmouth College/Exome-mutation-calls/", pattern = "\\.csv$")                          # create chr vector of all file names with .txt extension

#     1. Run 'isolate_snv_cols' fxn on 'filtered_variant_calls' files
#     2. Convert 'hgvs_c' and 'hgvs_p' to characters vectors
#     3. Write filtered .csv files

filtered_variant_calls
base_49_filt <- isolate_snv_cols(filtered_variant_calls[1])  

base_49_filt$hgvs_c <- as.character(base_49_filt$hgvs_c)   

base_50_filt <- isolate_snv_cols(filtered_variant_calls[2])
base_50_filt$hgvs_c <- as.character(base_50_filt$hgvs_c)

base_51_filt <- isolate_snv_cols(filtered_variant_calls[3])
base_51_filt$hgvs_c <- as.character(base_51_filt$hgvs_c)

recur_293_filt <- isolate_snv_cols(filtered_variant_calls[4])
recur_293_filt$hgvs_c <- as.character(recur_293_filt$hgvs_c)

recur_294_filt <- isolate_snv_cols(filtered_variant_calls[5])
recur_294_filt$hgvs_c <- as.character(recur_294_filt$hgvs_c)

recur_346_filt <- isolate_snv_cols(filtered_variant_calls[6])
recur_346_filt$hgvs_c <- as.character(recur_346_filt$hgvs_c)

recur_347_filt <- isolate_snv_cols(filtered_variant_calls[7])
recur_347_filt$hgvs_c <- as.character(recur_347_filt$hgvs_c)

recur_348_filt <- isolate_snv_cols(filtered_variant_calls[8])
recur_348_filt$hgvs_c <- as.character(recur_348_filt$hgvs_c)


# Merge & aggregate SNV data ----
base_merge = merge(base_49$hgvs_c, base_50$hgvs_c, all=TRUE) 
base_merge_2 = merge(base_merge_1, base_51, all=TRUE)
base_merge_2$counts <- 1                                                           # add column labeled 'counts'
base_snv_counts <- aggregate(base_merge_2$counts,list(base_merge_2$hgvs_c),sum)    # aggregate duplicates
colnames(base_snv_counts) <- c("hgvs_c", "counts")                                 # rename columns
write.csv(base_snv_counts, file = "Output_files/Baseline_SNV_counts.csv")                       # write file containing SNV counts

recur_merge_1 = merge(recur_293, recur_294, all=TRUE)
recur_merge_2 = merge(recur_merge_1, recur_346, all=TRUE)
recur_merge_3 = merge(recur_merge_2, recur_347, all=TRUE)
recur_merge_4 = merge(recur_merge_3, recur_348, all=TRUE)
recur_merge_4$counts <- 1                                                          # add column labeled 'counts'
recur_snv_counts <- aggregate(recur_merge_4$counts,list(recur_merge_4$hgvs_c),sum) # aggregate duplicates
colnames(recur_snv_counts) <- c("hgvs_c", "counts")                                # rename columns
write.csv(recur_snv_counts, file = "Output_files/Reccurence_SNV_counts.csv")                    # write file containing SNV counts


# Plot venn diagram ----
set1<-base_snv_counts$hgvs_c
set2<-recur_snv_counts$hgvs_c

venn.diagram(
  x = list(set1, set2),
  category.names = c("Baseline", "Recurrent"),
  filename = 'output_files/TIFFs/Venn_diagram/Venn.tiff',
  imagetype="tiff" ,
  resolution = 110,
  compression = "lzw",
  cex = 8,
  cat.cex=8,
  fill = c("#1B9E77", "#E7298A"),
  cat.pos = c(-1,2), 
  cat.default.pos = "text",
  cat.dist = c(0.16, 0.18), width = 3000
)


# SNv list unique in reccurrent tumors ----
unique_SNVs_reccurent <- anti_join (recur_snv_counts, base_snv_counts, by='hgvs_c')  
write.csv(unique_SNVs_reccurent, file = "Output_files/Unique_SNVs_in_reccurent_tumors.csv")



# -------------------------------------------------------------------------------
#
#    What biological significance does the SNV count data hold, if any?
#          -- technical 
#
# 

filter_fxn <- function(file_name) {
  new_file <- read.delim(file_name, na.strings = "")
  new_file[, 'hgvs_c'] <- as.character(new_file[, 'hgvs_c'])
  new_file[, 'hgvs_p'] <- as.character(new_file[, 'hgvs_p'])
  new_file_clean_filt <- new_file[na.omit(new_file$hgvs_p),]
}


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














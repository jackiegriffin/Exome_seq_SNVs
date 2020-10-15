# Load libraries ----
library(tidyverse)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)

# Upload data ----
variant_calls <- list.files(path = "C:/Users/jackie/OneDrive - Dartmouth College/Exome-mutation-calls/Input_variant_call_files/", pattern = "\\.txt$")                          # create chr vector of all file names with .txt extension
variant_calls <- paste('Input_variant_call_files/', variant_calls, sep = '')                  # paste folder name (containing .txt files) to reflect PATH from wd

# Write file process function ----
file_process_fxn <- function(file_name) {
  new_file <- read.delim(file_name)
  snv <- new_file[c('hgvs_c')]                                                     #'hgvs_c' = column containing SNV data 
}


# Run function ----
base_49 <- file_process_fxn(variant_calls[1])                                      # run
base_49$hgvs_c <- as.character(base_49$hgvs_c)                                     # convert data to character string
base_50 <- file_process_fxn(variant_calls[2])
base_50$hgvs_c <- as.character(base_50$hgvs_c)
base_51 <- file_process_fxn(variant_calls[3])
base_51$hgvs_c <- as.character(base_51$hgvs_c)

recur_293 <- file_process_fxn(variant_calls[4])
recur_293$hgvs_c <- as.character(recur_293$hgvs_c)
recur_294 <- file_process_fxn(variant_calls[5])
recur_294$hgvs_c <- as.character(recur_294$hgvs_c)
recur_346 <- file_process_fxn(variant_calls[6])
recur_346$hgvs_c <- as.character(recur_346$hgvs_c)
recur_347 <- file_process_fxn(variant_calls[7])
recur_347$hgvs_c <- as.character(recur_347$hgvs_c)
recur_348 <- file_process_fxn(variant_calls[8])
recur_348$hgvs_c <- as.character(recur_348$hgvs_c)


# Merge & aggregate SNV data ----
base_merge_1 = merge(base_49, base_50, all=TRUE) 
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

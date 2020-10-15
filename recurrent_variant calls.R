
# Load libraries ----
library(tidyverse)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)


# Functions ----
filter_out_synonymous <- function(file_name) {
  new_file <- read.delim(file_name, na.strings = "")
  new_file_filt <- na.omit(new_file, cols = "hgvs_p")
}

# Upload data ----
variant_calls <- list.files(path = "C:/Users/jackie/OneDrive - Dartmouth College/Exome-mutation-calls/Input_variant_call_files/", pattern = "\\.txt$")                          # create chr vector of all file names with .txt extension
variant_calls <- paste('Input_variant_call_files/', variant_calls, sep = '')  # paste folder name (containing .txt files) to reflect PATH from wd

#     1. Run 'filter_out_synonymous' fxn on 'variant_calls' files
#     2. Convert 'hgvs_c' and 'hgvs_p' to characters vectors
#     3. Write filtered .csv files

base_49 <- filter_out_synonymous(variant_calls[1])
base_49$hgvs_c <- as.character(base_49$hgvs_c)                                     
base_49$hgvs_p <- as.character(base_49$hgvs_p)       

base_50 <- filter_out_synonymous(variant_calls[2])
base_50$hgvs_c <- as.character(base_50$hgvs_c)                                     
base_50$hgvs_p <- as.character(base_50$hgvs_p)     

base_51 <- filter_out_synonymous(variant_calls[3])
base_51$hgvs_c <- as.character(base_51$hgvs_c)                                     
base_51$hgvs_p <- as.character(base_51$hgvs_p)          

recur_293 <- filter_out_synonymous(variant_calls[4])
recur_293$hgvs_c <- as.character(recur_293$hgvs_c)                                     
recur_293$hgvs_p <- as.character(recur_293$hgvs_p)      

recur_294 <- filter_out_synonymous(variant_calls[5])
recur_294$hgvs_c <- as.character(recur_294$hgvs_c)                                     
recur_294$hgvs_p <- as.character(recur_294$hgvs_p)     

recur_346 <- filter_out_synonymous(variant_calls[6])
recur_346$hgvs_c <- as.character(recur_346$hgvs_c)                                     
recur_346$hgvs_p <- as.character(recur_346$hgvs_p)       

recur_347 <- filter_out_synonymous(variant_calls[7])
recur_347$hgvs_c <- as.character(recur_347$hgvs_c)                                     
recur_347$hgvs_p <- as.character(recur_347$hgvs_p)       

recur_348 <- filter_out_synonymous(variant_calls[8])
recur_348$hgvs_c <- as.character(recur_348$hgvs_c)                                     
recur_348$hgvs_p <- as.character(recur_348$hgvs_p)      


# Subset mutations present in all bio replicates
genes_in_all_3_base <- Reduce(intersect, list(base_49$hgvs_c,base_50$hgvs_c,base_51$hgvs_c))
genes_in_all_5_recur <- Reduce(intersect, list(recur_293$hgvs_c,
                                                     recur_294$hgvs_c,
                                                     recur_346$hgvs_c,
                                                     recur_347$hgvs_c,
                                                     recur_348$hgvs_c))
#base_49_short <- base_49[base_49$hgvs_c %in% genes_in_all_3, ]


# Venn diagram ----
set1<-genes_in_all_3_base
set2<-genes_in_all_5_recur

venn.diagram(
  x = list(set1, set2),
  category.names = c("Baseline", "Recurrent"),
  filename = 'Venn.tiff',
  imagetype="tiff" ,
  resolution = 110,
  compression = "lzw",
  cex = 8,
  cat.cex=8,
  fill = c("#1B9E77", "#E7298A"),
  cat.pos = c(-3,2), 
  cat.default.pos = "text",
  cat.dist = c(0.16, 0.17), width = 3000
)


# Subset unique snvs in reccurent tumors----
genes_in_all_3_base <- as.data.frame(genes_in_all_3_base)
colnames(genes_in_all_3_base) <- c("hgvs_c")
genes_in_all_5_recur <- as.data.frame(genes_in_all_5_recur)
colnames(genes_in_all_5_recur) <- c("hgvs_c")
SNVs_unique_to_reccurent <- anti_join (genes_in_all_5_recur, genes_in_all_3_base, by='hgvs_c')  
write.csv(SNVs_unique_to_reccurent, file = "SNVs_unique_to_reccurent.csv")







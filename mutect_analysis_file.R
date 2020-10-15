library(dplyr)  


## have not in and mutect_p

twm_17_348 <- read.delim('upload_files/TWM_17_348_mutect_b37_ann.tsv',
                         stringsAsFactors = FALSE, 
                         sep = '\t', header = TRUE) # upload .ann.tsv file from mutect, snpeff and vcf2tsv
twm_17_348 <- mutect_process(twm_17_348) # run function




### stop before filtering and see if we find aa change in PIK3CA missense variant ----

twm_17_348 <-filter(twm_17_348, FILTER == "PASS") # only keep rows with filter column = pass
twm_17_348<-twm_17_348[!(twm_17_348$effect=="synonymous_variant"),] # remove synonymous variants


# write file to send to todd ---
write.table(twm_17_348, file ="twm_17_348.txt", sep = "\t", row.names = FALSE)

# ---------------------------------------------------------------------------------------------------------------

twm_17_348 <- read.delim('todd mutect files/twm_17_049_variants_filtered.txt',
                         stringsAsFactors = FALSE, 
                         sep = '\t', header = TRUE) # upload .ann.tsv file from mutect, snpeff and vcf2tsv
## PIK3CA variant encodes an amino acid change from 


twm_17_348 <- read.delim('',
                         stringsAsFactors = FALSE, 
                         sep = '\t', header = TRUE)

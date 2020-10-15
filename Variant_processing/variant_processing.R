library(dplyr)

#################################################################################
#                                                                               #
#    Somatic variants were called using gatk-mutect2-snpeff-vcf2tsv pipeline    #
#                                                                               #
#              BASH SCRIPT = 'RNA_Exomes.sh'                                    #
#              VARIANT PROCESSING = 'mutect_process' function                   #
#                                                                               #
#################################################################################


twm_17_049 <- read.delim('upload_files/TWM_17_049_mutect_b37_ann.tsv',
                         stringsAsFactors = FALSE, sep = '\t', header = TRUE) # load .ann.tsv files from gatk-mutect2 output
twm_17_049 <- mutect_process(twm_17_049) # run fxn
twm_17_049 <-filter(twm_17_049, FILTER == "PASS") # only keep rows with filter column = pass
twm_17_049 <-twm_17_049[!(twm_17_049$effect=="synonymous_variant"),] # remove synonymous variants
write.table(twm_17_049, file ="twm_17_049.txt", sep = "\t", row.names = FALSE) # write file to send to todd ---



## simple exclusion function ----
'%!in%' <- function(x,y)!('%in%'(x,y))

## function for processing mutect calls ----
mutect_process <- function(mutect_calls, sample_type = 'tumor') {
  
  mutect_info <- mutect_calls$ANN # just the annotation info from snpEff. 'ANN' stands for annotations
  info_split <- strsplit(mutect_info, '\\|') # split them by | and now a list
  
  mutant_allele <- rep(NA, length(info_split))
  for (i in 1:length(mutant_allele)) {
    mutant_allele[i] <- info_split[[i]][1]
  }
  mutect_calls$mutant_allele <- mutant_allele
  
  effect <- rep(NA, length(info_split))
  for (i in 1:length(effect)) {
    effect[i] <- info_split[[i]][2]
  }
  mutect_calls$effect <- effect
  
  impact <- rep(NA, length(info_split))
  for (i in 1:length(impact)) {
    impact[i] <- info_split[[i]][3]
  }
  mutect_calls$impact <- impact
  
  gene_name <- rep(NA, length(info_split))
  for (i in 1:length(gene_name)) {
    gene_name[i] <- info_split[[i]][4]
  }
  mutect_calls$gene_name <- gene_name
  
  gene_id <- rep(NA, length(info_split))
  for (i in 1:length(gene_id)) {
    gene_id[i] <- info_split[[i]][5]
  }
  mutect_calls$gene_id <- gene_id
  
  feature_type <- rep(NA, length(info_split))
  for (i in 1:length(feature_type)) {
    feature_type[i] <- info_split[[i]][6]
  }
  mutect_calls$feature_type <- feature_type
  
  feature_id <- rep(NA, length(info_split))
  for (i in 1:length(feature_id)) {
    feature_id[i] <- info_split[[i]][7]
  }
  mutect_calls$feature_id <- feature_id
  
  transcript_biotype <- rep(NA, length(info_split))
  for (i in 1:length(transcript_biotype)) {
    transcript_biotype[i] <- info_split[[i]][8]
  }
  mutect_calls$transcript_biotype <- transcript_biotype
  
  rank_total <- rep(NA, length(info_split))
  for (i in 1:length(rank_total)) {
    rank_total[i] <- info_split[[i]][9]
  }
  mutect_calls$rank_total <- rank_total
  
  hgvs_c <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_c)) {
    hgvs_c[i] <- info_split[[i]][10]
  }
  mutect_calls$hgvs_c <- hgvs_c
  
  hgvs_p <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_p)) {
    hgvs_p[i] <- info_split[[i]][11]
  }
  mutect_calls$hgvs_p <- hgvs_p
  
  cdna_pos <- rep(NA, length(info_split))
  for (i in 1:length(cdna_pos)) {
    cdna_pos[i] <- info_split[[i]][12]
  }
  mutect_calls$cdna_pos <- cdna_pos
  
  cds_pos_cds_len <- rep(NA, length(info_split))
  for (i in 1:length(cds_pos_cds_len)) {
    cds_pos_cds_len[i] <- info_split[[i]][13]
  }
  mutect_calls$cds_pos_cds_len <- cds_pos_cds_len
  
  protein_pos <- rep(NA, length(info_split))
  for (i in 1:length(protein_pos)) {
    protein_pos[i] <- info_split[[i]][14]
  }
  mutect_calls$protein_pos <- protein_pos
  
  dist_to_feat <- rep(NA, length(info_split))
  for (i in 1:length(dist_to_feat)) {
    dist_to_feat[i] <- info_split[[i]][15]
  }
  mutect_calls$dist_to_feat <- dist_to_feat
  
  #all calls with trusight genes and add location and remove ANN column
  mutect_calls$X.CHROM <- gsub('chr', '', mutect_calls$X.CHROM)
  mutect_calls$location <- paste0(mutect_calls$X.CHROM, ':', mutect_calls$POS)
  mutect_calls <- mutect_calls[, -which(names(mutect_calls) %in% c('ANN', 'QUAL', 'TLOD', 'gene_id', 'feature_id', 'cdna_pos'))]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intron_variant', ]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intergenic_region', ]
  mutect_calls$AF <- as.numeric(mutect_calls$AF)
  mutect_calls <- mutect_calls[!is.na(mutect_calls$AF), ]
  mutect_calls <- mutect_calls[nchar(mutect_calls$REF) == 1, ]
  mutect_calls <- mutect_calls[nchar(mutect_calls$ALT) == 1, ]
  mutect_calls <- mutect_calls[mutect_calls$CONTQ >= 20, ]
  #mutect_calls <- mutect_calls[mutect_calls$ECNT == 1, ]
  mutect_calls <- mutect_calls[mutect_calls$GERMQ >= 20, ]
  mutect_calls <- mutect_calls[mutect_calls$SEQQ >= 20, ]
  #if (sample_type == 'plasma') {
   # mutect_calls <- mutect_calls[mutect_calls$AF >= 0.05, ]
  #}
  #else {
    mutect_calls <- mutect_calls[mutect_calls$AF >= 0.10, ]
  #}
  mutect_calls <- mutect_calls[, which(colnames(mutect_calls) %in% c('location', 'hgvs_c', 'hgvs_p', 'gene_name', 'impact', 'effect', 'AF', 'DP', 'DP.1', 'FILTER'))]
  return(mutect_calls)
}

write.table(twm_17_049, file ="twm_17_049.txt", sep = "\t", row.names = FALSE)


## function to draw tile figure ----
draw_tiles <- function(sample_1, sample_2, sample_3, sample_4 = FALSE, plasma_sample, row_labels = NULL, fill_NA = FALSE, 
                       color_scheme = c('green/purple', 'red/blue'), col_labels = c('LOCATION', 'GENE')) { 
  patient_met_pool <- unique(c(sample_1$location, sample_2$location, sample_3$location)) #get all of the locations found in the tumor samples
  
  if (!missing(sample_4)) { #add in sample_4 if given
    patient_met_pool <- unique(c(patient_met_pool, sample_4$location))
  }
  
  plasma_found <- plasma_sample[plasma_sample$location %in% patient_met_pool, ] #finds locations also found in plasma
  plasma_found_vars <- (plasma_found$location) #these locations to be seen in tumor samples
  
  # subset to variants found in plasma and take a look
  sample_1_pooled <- sample_1[sample_1$location %in% plasma_found_vars, ]
  
  sample_2_pooled <- sample_2[sample_2$location %in% plasma_found_vars, ]
  
  sample_3_pooled <- sample_3[sample_3$location %in% plasma_found_vars, ]
  
  if (!missing(sample_4)) { #add sample_4 if given
    sample_4_pooled <- sample_4[sample_4$location %in% plasma_found_vars, ]
  }
  
  plasma_sample_pooled <- plasma_sample[plasma_sample$location %in% plasma_found_vars, ]
  
  # subset these to just tumor AF for tile figure
  sample_1_pooled <- sample_1_pooled[, c('location', 'AF')]
  colnames(sample_1_pooled) <- c('location', 'AF_sample_1')

  sample_2_pooled <- sample_2_pooled[, c('location', 'AF')]
  colnames(sample_2_pooled) <- c('location', 'AF_sample_2')

  sample_3_pooled <- sample_3_pooled[, c('location', 'AF')]
  colnames(sample_3_pooled) <- c('location', 'AF_sample_3')
  
  if (!missing(sample_4)) {
    sample_4_pooled <- sample_4_pooled[, c('location', 'AF')]
    colnames(sample_4_pooled) <- c('location', 'AF_sample_4')
  }

  plasma_sample_pooled <- plasma_sample_pooled[, c('location', 'AF')]
  colnames(plasma_sample_pooled) <- c('location', 'AF_plasma')
  
  # put them together in one dataframe
  muts_pooled <- merge(sample_1_pooled, sample_2_pooled, by = 'location', all = TRUE)
  muts_pooled <- merge(muts_pooled, sample_3_pooled, by = 'location', all = TRUE)
  
  if (!missing(sample_4)) {
    muts_pooled <- merge(muts_pooled, sample_4_pooled, by = 'location', all = TRUE)
  }
  
  muts_pooled <- merge(muts_pooled, plasma_sample_pooled, by = 'location', all = TRUE)
  muts_pooled <- muts_pooled[order(muts_pooled$AF_plasma, decreasing = TRUE), ]
  
  rownames(muts_pooled) <- muts_pooled$location #make location the row names
  row_muts <- rownames(muts_pooled)
  muts_pooled <- muts_pooled[, -1] #remove location variable from df
  
  muts_pooled <- as.data.frame(muts_pooled)
  rownames(muts_pooled) <- row_muts
  
  #plot 
  
  tumors <- c(muts_pooled$AF_sample_1, muts_pooled$AF_sample_2, muts_pooled$AF_sample_3) #everything but plasma
  
  if (!missing(sample_4)) {
    tumors <- c(tumors, muts_pooled$AF_sample_4)
  }
  
  #set colors, plasma has its own reds, pooled tumors blues
  cell_cols<-rep("#000000",dim(muts_pooled)[1] * dim(muts_pooled)[2]) #create matrix to hold colors for graphing
  
    
  if (color_scheme == 'red/blue') {
    # plasma colors
    if (!missing(sample_4)) {
      cell_cols[(length(plasma_found_vars)*4) + 1:(length(plasma_found_vars)*5)] <- color.scale(muts_pooled[, ncol(muts_pooled)], extremes = c('maroon1', 'maroon4'), na.color = '#ffffff')
    }
    else {
      cell_cols[(length(plasma_found_vars)*3) + 1:(length(plasma_found_vars)*4)] <- color.scale(muts_pooled[, ncol(muts_pooled)], extremes = c('maroon1', 'maroon4'), na.color = '#ffffff')
    }
  
    # tumor colors
    if (!missing(sample_4)) {
      cell_cols[1:(length(plasma_found_vars)*4)] <- color.scale(tumors, extremes = c('turquoise1', 'turquoise4'), na.color = '#ffffff')
    }
    else {
      cell_cols[1:(length(plasma_found_vars)*3)] <- color.scale(tumors, extremes = c('turquoise1', 'turquoise4'), na.color = '#ffffff')
    }
  }
  
  else if (color_scheme == 'green/purple') {
    # plasma colors
    if (!missing(sample_4)) {
      cell_cols[(length(plasma_found_vars)*4) + 1:(length(plasma_found_vars)*5)] <- color.scale(muts_pooled[, ncol(muts_pooled)], extremes = c('seagreen1', 'seagreen4'), na.color = '#ffffff')
    }
    else {
      cell_cols[(length(plasma_found_vars)*3) + 1:(length(plasma_found_vars)*4)] <- color.scale(muts_pooled[, ncol(muts_pooled)], extremes = c('seagreen1', 'seagreen4'), na.color = '#ffffff')
    }
    
    # tumor colors
    if (!missing(sample_4)) {
      cell_cols[1:(length(plasma_found_vars)*4)] <- color.scale(tumors, extremes = c('darkorchid1', 'darkorchid4'), na.color = '#ffffff')
    }
    else {
      cell_cols[1:(length(plasma_found_vars)*3)] <- color.scale(tumors, extremes = c('darkorchid1', 'darkorchid4'), na.color = '#ffffff')
    }
  }
  
  cell_cols <- matrix(cell_cols, nrow = length(plasma_found_vars), byrow = FALSE) #put into matrix form
  
  pooled_t <- data.frame(t(muts_pooled))
  pooled_t <- pooled_t[c(ncol(muts_pooled), 1:(ncol(muts_pooled) - 1)), ] 
  
  cell_cols <- t(cell_cols)
  cell_cols <- cell_cols[c(ncol(muts_pooled), 1:(ncol(muts_pooled) - 1)), ] 
  
  # plot it
  # extra space
  par(mar=c(6,5.5,6,2.1))
  #par(mar=c(6,15.5,6,12.1))
  color2D.matplot(pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE) #SET AXIS LABELS AS OPTIONS AND ADD MAIN TITLE OPTION
  
  # add column labels SET OPTIONS FOR THESE (ALL, PROTEIN ONLY, LOCATION, GENE ONLY)!!!!!!!!!!!!!!!!!!!!
  if (col_labels == 'LOCATION') {
    axis(3, at = (1:ncol(pooled_t)) - 0.6, labels = rownames(muts_pooled), tick = FALSE, cex.axis = 0.7, las = 2, font = 2)
  }
  
  else if (col_labels == 'GENE') {
    axis(3, at = (1:ncol(pooled_t)) - 0.6, labels = muts_pooled$gene_name[muts_pooled$locations == rownames(muts_pooled)], tick = FALSE, cex.axis = 0.7, las = 2, font = 2)
  }
  
  # add row labels 
  if (!missing(sample_4)) {
    axis(2, at = c(0.5, 1.5, 2.5, 3.5, 4.5), labels = rev(row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)
  }
  else {
    axis(2, at = c(0.5, 1.5, 2.5, 3.5), labels = rev(row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)
  }
  
  # add NA markers
  if (fill_NA == TRUE) {
    if (!missing(sample_4)) {
    # sample_1 points
      points(x = which(is.na(muts_pooled$AF_sample_1)) - 0.5, 
            y = rep(3.5, sum(is.na(muts_pooled$AF_sample_1))), 
            pch = 16)
    # sample_2 points
      points(x = which(is.na(muts_pooled$AF_sample_2)) - 0.5, 
            y = rep(2.5, sum(is.na(muts_pooled$AF_sample_2))), 
            pch = 16)
    # sample_3 points
      points(x = which(is.na(muts_pooled$AF_sample_3)) - 0.5, 
            y = rep(1.5, sum(is.na(muts_pooled$AF_sample_3))), 
            pch = 16)
    # sample_4 points
      points(x = which(is.na(muts_pooled$AF_sample_4)) - 0.5, 
             y = rep(0.5, sum(is.na(muts_pooled$AF_sample_4))), 
             pch = 16)
    }
    else {
      # sample_1 points
      points(x = which(is.na(muts_pooled$AF_sample_1)) - 0.5, 
             y = rep(2.5, sum(is.na(muts_pooled$AF_sample_1))), 
             pch = 16)
      # sample_2 points
      points(x = which(is.na(muts_pooled$AF_sample_2)) - 0.5, 
             y = rep(1.5, sum(is.na(muts_pooled$AF_sample_2))), 
             pch = 16)
      # sample_3 points
      points(x = which(is.na(muts_pooled$AF_sample_3)) - 0.5, 
             y = rep(0.5, sum(is.na(muts_pooled$AF_sample_3))), 
             pch = 16)
    }
  }
}


mutect_process_indel <- function(mutect_calls, sample_type = 'tumor') {
  
  mutect_info <- mutect_calls$ANN # just the annotation info from snpEff
  info_split <- strsplit(mutect_info, '\\|') # split them by | and now a list
  
  mutant_allele <- rep(NA, length(info_split))
  for (i in 1:length(mutant_allele)) {
    mutant_allele[i] <- info_split[[i]][1]
  }
  mutect_calls$mutant_allele <- mutant_allele
  
  effect <- rep(NA, length(info_split))
  for (i in 1:length(effect)) {
    effect[i] <- info_split[[i]][2]
  }
  mutect_calls$effect <- effect
  
  impact <- rep(NA, length(info_split))
  for (i in 1:length(impact)) {
    impact[i] <- info_split[[i]][3]
  }
  mutect_calls$impact <- impact
  
  gene_name <- rep(NA, length(info_split))
  for (i in 1:length(gene_name)) {
    gene_name[i] <- info_split[[i]][4]
  }
  mutect_calls$gene_name <- gene_name
  
  gene_id <- rep(NA, length(info_split))
  for (i in 1:length(gene_id)) {
    gene_id[i] <- info_split[[i]][5]
  }
  mutect_calls$gene_id <- gene_id
  
  feature_type <- rep(NA, length(info_split))
  for (i in 1:length(feature_type)) {
    feature_type[i] <- info_split[[i]][6]
  }
  mutect_calls$feature_type <- feature_type
  
  feature_id <- rep(NA, length(info_split))
  for (i in 1:length(feature_id)) {
    feature_id[i] <- info_split[[i]][7]
  }
  mutect_calls$feature_id <- feature_id
  
  transcript_biotype <- rep(NA, length(info_split))
  for (i in 1:length(transcript_biotype)) {
    transcript_biotype[i] <- info_split[[i]][8]
  }
  mutect_calls$transcript_biotype <- transcript_biotype
  
  rank_total <- rep(NA, length(info_split))
  for (i in 1:length(rank_total)) {
    rank_total[i] <- info_split[[i]][9]
  }
  mutect_calls$rank_total <- rank_total
  
  hgvs_c <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_c)) {
    hgvs_c[i] <- info_split[[i]][10]
  }
  mutect_calls$hgvs_c <- hgvs_c
  
  hgvs_p <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_p)) {
    hgvs_p[i] <- info_split[[i]][11]
  }
  mutect_calls$hgvs_p <- hgvs_p
  
  cdna_pos <- rep(NA, length(info_split))
  for (i in 1:length(cdna_pos)) {
    cdna_pos[i] <- info_split[[i]][12]
  }
  mutect_calls$cdna_pos <- cdna_pos
  
  cds_pos_cds_len <- rep(NA, length(info_split))
  for (i in 1:length(cds_pos_cds_len)) {
    cds_pos_cds_len[i] <- info_split[[i]][13]
  }
  mutect_calls$cds_pos_cds_len <- cds_pos_cds_len
  
  protein_pos <- rep(NA, length(info_split))
  for (i in 1:length(protein_pos)) {
    protein_pos[i] <- info_split[[i]][14]
  }
  mutect_calls$protein_pos <- protein_pos
  
  dist_to_feat <- rep(NA, length(info_split))
  for (i in 1:length(dist_to_feat)) {
    dist_to_feat[i] <- info_split[[i]][15]
  }
  mutect_calls$dist_to_feat <- dist_to_feat
  
  #all calls with trusight genes and add location and remove ANN column
  mutect_calls$X.CHROM <- gsub('chr', '', mutect_calls$X.CHROM)
  mutect_calls$location <- paste0(mutect_calls$X.CHROM, ':', mutect_calls$POS)
  mutect_calls <- mutect_calls[, -which(names(mutect_calls) %in% c('ANN', 'QUAL', 'FILTER', 'TLOD', 'gene_id', 'feature_id', 'cdna_pos'))]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intron_variant', ]
  mutect_calls <- mutect_calls[mutect_calls$effect != 'intergenic_region', ]
  mutect_calls$AF <- as.numeric(mutect_calls$AF)
  mutect_calls <- mutect_calls[!is.na(mutect_calls$AF), ]
  #mutect_calls <- mutect_calls[nchar(mutect_calls$REF) == 1, ]
  #mutect_calls <- mutect_calls[nchar(mutect_calls$ALT) == 1, ]
  mutect_calls <- mutect_calls[mutect_calls$CONTQ >= 20, ]
  #mutect_calls <- mutect_calls[mutect_calls$ECNT == 1, ]
  mutect_calls <- mutect_calls[mutect_calls$GERMQ >= 20, ]
  mutect_calls <- mutect_calls[mutect_calls$SEQQ >= 20, ]
  if (sample_type == 'plasma') {
    mutect_calls <- mutect_calls[mutect_calls$AF >= 0.01, ]
  }
  else {
    mutect_calls <- mutect_calls[mutect_calls$AF >= 0.05, ]
  }
  return(mutect_calls)
}

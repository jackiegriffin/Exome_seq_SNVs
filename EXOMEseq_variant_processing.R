# load libraries ----
      library(dplyr)

################################################################################
#                                                                              #
#    PLATFORM : EXOME-seq                                                      #
#    PIPELINE : SNV calls using gatk/mutect2-snpeff-vcf2tsv ('EXOMEseq.sh')    #
#    INPUT (pipeline output) : mutect_filtered_ann.tsv files                   #
#    Variant processing using custom R function ('mutect_process.R')           #
#    OUTPUT* : somatic missense mutations unique to recurrent breast tumors    #
#         *Todd wants RAW and FILTERED hgvs_c missense mutation LIST.txt       #
#         * See NOTES TO SELF at end of script                                 #
#                                                                              #
################################################################################

# ************* functions  ************************************************************----

  # Simple exclusion function ----
        '%!in%' <- function(x,y)!('%in%'(x,y))
        
  # Function for processing mutect calls ----
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
        
  # Function for filtering out synonymous variants ----
        filter_out_synonymous <- function(file_name) {
          new_file <- read.delim(file_name, na.strings = "")
          new_file_filt <- na.omit(new_file, cols = "hgvs_p")
          new_fill_filt_pass <- filter(new_file_filt, FILTER== "PASS")
          new_fill_filt_pass_syn <-new_fill_filt_pass[!(new_fill_filt_pass$effect=="synonymous_variant"),] # remove synonymous variants
          new_fil_filt_allele <- new_fill_filt_pass_syn[(new_fill_filt_pass_syn$AF >=0.1 & new_fill_filt_pass_syn$AF <=0.5),]
          new_fil_filt_allele_impact <- filter(new_fil_filt_allele, impact=="HIGH"| impact=="MODERATE")
        }
        
# ****************************************************************************************----
    
# Load files and run fxn ----
      
        tsv_files <- list.files(path="mutect_ann.tsv_files",  # list .tsv files
                     pattern = "*.tsv$",
                     full.names = TRUE)  
      
        lst <- vector("list", length(tsv_files))  # create empty list
      
        for(i in 1:length(tsv_files)) { # read .tsv files into list
          lst[[i]] <- read.delim(tsv_files[i], 
                      sep ='\t', header =TRUE, 
                      stringsAsFactors = FALSE)
        }
        
        names(lst) <- c(tsv_files) # re-name list as file names
        lst_mutect_processed <- lapply(lst, mutect_process) # run fxn
        df <- lapply(lst_mutect_processed, as.data.frame) # covert to df
        
        
# Write variant calls for individual tumor samples (RAW OUTPUT = SEND TO TODD) ----
        
    # Control (primary)
        control_049 <- df$`mutect_ann.tsv_files/TWM_17_049_mutect_b37_ann.tsv`
          write.table(control_049, file ="Filtered_variant_files/mutect_processed/control_049_raw.txt", sep = "\t", row.names = FALSE)
        control_050 <- df$`mutect_ann.tsv_files/TWM_17_050_mutect_b37_ann.tsv`
          write.table(control_050, file ="Filtered_variant_files/mutect_processed/control_050_raw.txt", sep = "\t", row.names = FALSE)
        control_051 <- df$`mutect_ann.tsv_files/TWM_17_051_mutect_b37_ann.tsv`
          write.table(control_051, file ="Filtered_variant_files/mutect_processed/control_051_raw.txt", sep = "\t", row.names = FALSE)
    # Resistant (recurrent)  
        recurrent_238 <- df$`mutect_ann.tsv_files/TWM_17_238_mutect_b37_ann.tsv`
          write.table(recurrent_238, file ="Filtered_variant_files/mutect_processed/recurrent_238_raw.txt", sep = "\t", row.names = FALSE)
        recurrent_240 <- df$`mutect_ann.tsv_files/TWM_17_240_mutect_b37_ann.tsv`
          write.table(recurrent_240, file ="Filtered_variant_files/mutect_processed/recurrent_240_raw.txt", sep = "\t", row.names = FALSE)
        recurrent_346 <- df$`mutect_ann.tsv_files/TWM_17_346_mutect_b37_ann.tsv`
          write.table(recurrent_346, file ="Filtered_variant_files/mutect_processed/recurrent_346_raw.txt", sep = "\t", row.names = FALSE)
        recurrent_347 <- df$`mutect_ann.tsv_files/TWM_17_347_mutect_b37_ann.tsv`
          write.table(recurrent_347, file ="Filtered_variant_files/mutect_processed/recurrent_347_raw.txt", sep = "\t", row.names = FALSE)
        recurrent_348 <- df$`mutect_ann.tsv_files/TWM_17_348_mutect_b37_ann.tsv`
          write.table(recurrent_348, file ="Filtered_variant_files/mutect_processed/recurrent_348_raw.txt", sep = "\t", row.names = FALSE)
        

# Filter RAW data for non-synonymous and allele frequency ----  
    # Upload data 
        variant_calls_RAW <- list.files(path="mutect_processed/", pattern = "\\.txt$")                          
        variant_calls_RAW <- paste('mutect_processed/', variant_calls_RAW, sep = '')
        variant_calls_RAW
        
    # Control (primary)
        control_049_filt <- filter_out_synonymous(variant_calls_RAW[1])
        control_050_filt <- filter_out_synonymous(variant_calls_RAW[2])
        control_051_filt <- filter_out_synonymous(variant_calls_RAW[3])
    # Resistant (recurrent)  
        recurrent_238_filt <- filter_out_synonymous(variant_calls_RAW[4])
        recurrent_240_filt <- filter_out_synonymous(variant_calls_RAW[5])
        recurrent_346_filt <- filter_out_synonymous(variant_calls_RAW[6])
        recurrent_347_filt <- filter_out_synonymous(variant_calls_RAW[7])
        recurrent_348_filt <- filter_out_synonymous(variant_calls_RAW[8])

        
        
# Subset based on common variants across all biological replicates ----
        control_inters_genes <- as.data.frame(Reduce(intersect, list(control_049_filt$gene_name,
                                                                control_050_filt$gene_name,
                                                                control_051_filt$gene_name)))
        colnames(control_inters_genes) <- c("gene_name")
        
        recurrrent_inters_genes <- as.data.frame(Reduce(intersect, list(recurrent_238_filt$gene_name,
                                                                recurrent_240_filt$gene_name,
                                                                recurrent_346_filt$gene_name,
                                                                recurrent_347_filt$gene_name,
                                                                recurrent_348_filt$gene_name)))
        colnames(recurrrent_inters_genes) <- c("gene_name")
        
# Identify variants uniquely occuring in recurrent tumors ----
        recurrent_inters_UNIQUE <- anti_join (recurrrent_inters_genes, control_inters_genes, by='gene_name')  
        
# Subset based on unique variants in recurrent tumors only (FILTERED, UNIQUE OUTPUT = SEND TO TODD) ----        
        recurrent_238_filt_UNIQUE <- merge(recurrent_238_filt, recurrent_inters_UNIQUE['gene_name'])
          write.table(recurrent_238_filt_UNIQUE, file ="mutect_processed/filtered/unique/recurrent_238_filt_UNIQUE.txt", sep = "\t", row.names = FALSE)
        recurrent_240_filt_UNIQUE <- merge(recurrent_240_filt, recurrent_inters_UNIQUE['gene_name'])
          write.table(recurrent_240_filt_UNIQUE, file ="mutect_processed/filtered/unique/recurrent_240_filt_UNIQUE.txt", sep = "\t", row.names = FALSE)
        recurrent_346_filt_UNIQUE <- merge(recurrent_346_filt, recurrent_inters_UNIQUE['gene_name'])
          write.table(recurrent_346_filt_UNIQUE, file ="mutect_processed/filtered/unique/recurrent_346_filt_UNIQUE.txt", sep = "\t", row.names = FALSE)
        recurrent_347_filt_UNIQUE <- merge(recurrent_347_filt, recurrent_inters_UNIQUE['gene_name'])
          write.table(recurrent_347_filt_UNIQUE, file ="mutect_processed/filtered/unique/recurrent_347_filt_UNIQUE.txt", sep = "\t", row.names = FALSE)
        recurrent_348_filt_UNIQUE <- merge(recurrent_348_filt, recurrent_inters_UNIQUE['gene_name'])
          write.table(recurrent_348_filt_UNIQUE, file ="mutect_processed/filtered/unique/recurrent_348_filt_UNIQUE.txt", sep = "\t", row.names = FALSE)
          
# NOTES TO SELF ---
      # Same Exome-seq files were analyzed using:
          # CNVkit
          # Phylogenetic analysis
          # SigProfilerExtractor/plotter ( python3 env)
          # IGV -- original vcf2maf in DBI for mySQL clinical relevance query
          
        
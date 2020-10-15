# clean mutation analyses


## function for trimming mutect calls to trusight genes ----

trusight_from_mutect <- function(mutect_short) {
  ann_idx <- mutect_short$ANN == '' # index of mutations without annotations
  mutect_shorter <- mutect_short[!ann_idx, ] # remove any without annotation
  
  mutect_info <- mutect_shorter$ANN # just the annotation info from snpEff
  info_split <- strsplit(mutect_info, '\\|') # split them by | and now a list
  
  mutant_allele <- rep(NA, length(info_split))
  for (i in 1:length(mutant_allele)) {
    mutant_allele[i] <- info_split[[i]][1]
  }
  mutect_shorter$mutant_allele <- mutant_allele
  
  effect <- rep(NA, length(info_split))
  for (i in 1:length(effect)) {
    effect[i] <- info_split[[i]][2]
  }
  mutect_shorter$effect <- effect
  
  impact <- rep(NA, length(info_split))
  for (i in 1:length(impact)) {
    impact[i] <- info_split[[i]][3]
  }
  mutect_shorter$impact <- impact
  
  gene_name <- rep(NA, length(info_split))
  for (i in 1:length(gene_name)) {
    gene_name[i] <- info_split[[i]][4]
  }
  mutect_shorter$gene_name <- gene_name
  
  gene_id <- rep(NA, length(info_split))
  for (i in 1:length(gene_id)) {
    gene_id[i] <- info_split[[i]][5]
  }
  mutect_shorter$gene_id <- gene_id
  
  feature_type <- rep(NA, length(info_split))
  for (i in 1:length(feature_type)) {
    feature_type[i] <- info_split[[i]][6]
  }
  mutect_shorter$feature_type <- feature_type
  
  feature_id <- rep(NA, length(info_split))
  for (i in 1:length(feature_id)) {
    feature_id[i] <- info_split[[i]][7]
  }
  mutect_shorter$feature_id <- feature_id
  
  transcript_biotype <- rep(NA, length(info_split))
  for (i in 1:length(transcript_biotype)) {
    transcript_biotype[i] <- info_split[[i]][8]
  }
  mutect_shorter$transcript_biotype <- transcript_biotype
  
  rank_total <- rep(NA, length(info_split))
  for (i in 1:length(rank_total)) {
    rank_total[i] <- info_split[[i]][9]
  }
  mutect_shorter$rank_total <- rank_total
  
  hgvs_c <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_c)) {
    hgvs_c[i] <- info_split[[i]][10]
  }
  mutect_shorter$hgvs_c <- hgvs_c
  
  hgvs_p <- rep(NA, length(info_split))
  for (i in 1:length(hgvs_p)) {
    hgvs_p[i] <- info_split[[i]][11]
  }
  mutect_shorter$hgvs_p <- hgvs_p
  
  cdna_pos <- rep(NA, length(info_split))
  for (i in 1:length(cdna_pos)) {
    cdna_pos[i] <- info_split[[i]][12]
  }
  mutect_shorter$cdna_pos <- cdna_pos
  
  cds_pos_cds_len <- rep(NA, length(info_split))
  for (i in 1:length(cds_pos_cds_len)) {
    cds_pos_cds_len[i] <- info_split[[i]][13]
  }
  mutect_shorter$cds_pos_cds_len <- cds_pos_cds_len
  
  protein_pos <- rep(NA, length(info_split))
  for (i in 1:length(protein_pos)) {
    protein_pos[i] <- info_split[[i]][14]
  }
  mutect_shorter$protein_pos <- protein_pos
  
  dist_to_feat <- rep(NA, length(info_split))
  for (i in 1:length(dist_to_feat)) {
    dist_to_feat[i] <- info_split[[i]][15]
  }
  mutect_shorter$dist_to_feat <- dist_to_feat
  #how many genes in trusight panel?
  trusight <- intersect(mutect_shorter$gene_name, trusight_genes)
  
  #all calls with trusight genes
  mutect_trusight <- mutect_shorter[mutect_shorter$gene_name %in% trusight, ] 
  return(mutect_trusight)
}

## trusight 170 panel genes ----
trusight_genes <- c('AKT1', 'BRIP1', 'CREBBP', 'FANCI', 'FGFR2', 'JAK3', 'MSH3', 'PALB2', 'RAD51D', 'TSC1', 'AKT2', 'BTK', 'CSF1R', 'FANCL', 'FGFR3', 'KDR', 'MSH6', 'PDGFRA', 'RAD54L', 'TSC2',
                    'AKT3', 'CARD11', 'CTNNB1', 'FBXW7', 'FGFR4', 'KIT', 'MTOR', 'PDGFRB', 'RB1', 'VHL', 'ALK', 'CCND1', 'DDR2', 'FGF1', 'FLT1', 'KMT2A', 'MLL', 'MUTYH', 'PIK3CA', 'RET', 'XRCC2',
                    'APC', 'CCND2', 'DNMT3A', 'FGF2', 'FLT3', 'KRAS', 'MYC', 'PIK3CB', 'RICTOR', 'AR', 'CCNE1', 'EGFR', 'FGF3', 'FOXL2', 'MAP2K1', 'MYCL1', 'PIK3CD', 'ROS1',
                    'ARID1A', 'CD79A', 'EP300', 'FGF4', 'GEN1', 'MAP2K2', 'MYCN', 'PIK3CG', 'RPS6KB1', 'ATM', 'CD79B', 'ERBB2', 'FGF5', 'GNA11', 'MCL1', 'MYD88', 'PIK3R1', 'SLX4',
                    'ATR', 'CDH1', 'ERBB3', 'FGF6', 'GNAQ', 'MDM2', 'NBN', 'PMS2', 'SMAD4', 'BAP1', 'CDK12', 'ERBB4', 'FGF7', 'GNAS', 'MDM4', 'NF1', 'PPP2R2A', 'SMARCB1',
                    'BARD1', 'CDK4', 'ERCC1', 'FGF8', 'HNF1A', 'MET', 'NOTCH1', 'PTCH1', 'SMO', 'BCL2', 'CDK6', 'ERCC2', 'FGF9', 'HRAS', 'MLH1', 'NOTCH2', 'PTEN', 'SRC',
                    'BCL6', 'CDKN2A', 'ERG', 'FGF10', 'IDH1', 'MLLT3', 'NOTCH3', 'PTPN11', 'STK11', 'BRAF', 'CEBPA', 'ESR1', 'FGF14', 'IDH2', 'MPL', 'NPM1', 'RAD51', 'TERT',
                    'BRCA1', 'CHEK1', 'EZH2', 'FGF23', 'INPP4B', 'MRE11A', 'NRAS', 'RAD51B', 'TET2', 'BRCA2', 'CHEK2', 'FAM175A', 'FGFR1', 'JAK2', 'MSH2', 'NRG1', 'RAD51C', 'TP53')

## PATIENT 10 ----
## liver 1----
# import data
pat_10_liver_1_ug <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #12953
pat_10_liver_1_varscan <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #6056
pat_10_liver_1_freebayes <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #112330

pat_10_liver_1_mutect <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #188803

# dummy variable for matching
pat_10_liver_1_ug$location <- paste0(pat_10_liver_1_ug$CHROM, ':', pat_10_liver_1_ug$POS)
pat_10_liver_1_varscan$location <- paste0(pat_10_liver_1_varscan$CHROM, ':', pat_10_liver_1_varscan$POS)
pat_10_liver_1_freebayes$location <- paste0(pat_10_liver_1_freebayes$CHROM, ':', pat_10_liver_1_freebayes$POS)

pat_10_liver_1_mutect$location <- paste0(pat_10_liver_1_mutect$CHROM, ':', pat_10_liver_1_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_1_ug_fb <- intersect(pat_10_liver_1_ug$location, pat_10_liver_1_freebayes$location) #11581
pat_10_liver_1_ug_vs <- intersect(pat_10_liver_1_ug$location, pat_10_liver_1_varscan$location) #2733
pat_10_liver_1_vs_fb <- intersect(pat_10_liver_1_varscan$location, pat_10_liver_1_freebayes$location) #2673

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_1_mutect_ug_fb <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_ug_fb) #6526
pat_10_liver_1_mutect_ug_vs <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_ug_vs) #658
pat_10_liver_1_mutect_vs_fb <- intersect(pat_10_liver_1_mutect$location, pat_10_liver_1_vs_fb) #602

#combine them
pat_10_liver_1_all <- c(pat_10_liver_1_mutect_ug_fb, pat_10_liver_1_mutect_ug_vs, pat_10_liver_1_mutect_vs_fb)

#subset mutect calls
pat_10_liver_1_mutect_short <- pat_10_liver_1_mutect[pat_10_liver_1_mutect$location %in% pat_10_liver_1_all, ] #6680

#trim to trusight gene calls
pat_10_liver_1_mutect_trusight <- trusight_from_mutect(pat_10_liver_1_mutect_short) #17

## liver 2a----
# import data
pat_10_liver_2a_ug <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_2a_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #9193
pat_10_liver_2a_varscan <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_2a_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3947
pat_10_liver_2a_freebayes <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_2a_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #85821

pat_10_liver_2a_mutect <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_2a_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #124767

# dummy variable for matching
pat_10_liver_2a_ug$location <- paste0(pat_10_liver_2a_ug$CHROM, ':', pat_10_liver_2a_ug$POS)
pat_10_liver_2a_varscan$location <- paste0(pat_10_liver_2a_varscan$CHROM, ':', pat_10_liver_2a_varscan$POS)
pat_10_liver_2a_freebayes$location <- paste0(pat_10_liver_2a_freebayes$CHROM, ':', pat_10_liver_2a_freebayes$POS)

pat_10_liver_2a_mutect$location <- paste0(pat_10_liver_2a_mutect$CHROM, ':', pat_10_liver_2a_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_2a_ug_fb <- intersect(pat_10_liver_2a_ug$location, pat_10_liver_2a_freebayes$location) #8146
pat_10_liver_2a_ug_vs <- intersect(pat_10_liver_2a_ug$location, pat_10_liver_2a_varscan$location) #2182
pat_10_liver_2a_vs_fb <- intersect(pat_10_liver_2a_varscan$location, pat_10_liver_2a_freebayes$location) #2073

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_2a_mutect_ug_fb <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_ug_fb) #3639
pat_10_liver_2a_mutect_ug_vs <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_ug_vs) #393
pat_10_liver_2a_mutect_vs_fb <- intersect(pat_10_liver_2a_mutect$location, pat_10_liver_2a_vs_fb) #354

#combine them
pat_10_liver_2a_all <- c(pat_10_liver_2a_mutect_ug_fb, pat_10_liver_2a_mutect_ug_vs, pat_10_liver_2a_mutect_vs_fb)

#subset mutect calls
pat_10_liver_2a_mutect_short <- pat_10_liver_2a_mutect[pat_10_liver_2a_mutect$location %in% pat_10_liver_2a_all, ] #3728

#trim to trusight gene calls
pat_10_liver_2a_mutect_trusight <- trusight_from_mutect(pat_10_liver_2a_mutect_short) #10

## liver 5----
# import data
pat_10_liver_5_ug <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_5_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #12980
pat_10_liver_5_varscan <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_5_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4821
pat_10_liver_5_freebayes <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_5_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #118457

pat_10_liver_5_mutect <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_liver_5_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #193162

# dummy variable for matching
pat_10_liver_5_ug$location <- paste0(pat_10_liver_5_ug$CHROM, ':', pat_10_liver_5_ug$POS)
pat_10_liver_5_varscan$location <- paste0(pat_10_liver_5_varscan$CHROM, ':', pat_10_liver_5_varscan$POS)
pat_10_liver_5_freebayes$location <- paste0(pat_10_liver_5_freebayes$CHROM, ':', pat_10_liver_5_freebayes$POS)

pat_10_liver_5_mutect$location <- paste0(pat_10_liver_5_mutect$CHROM, ':', pat_10_liver_5_mutect$POS)

#pairwise intersections for the 3
pat_10_liver_5_ug_fb <- intersect(pat_10_liver_5_ug$location, pat_10_liver_5_freebayes$location) #11857
pat_10_liver_5_ug_vs <- intersect(pat_10_liver_5_ug$location, pat_10_liver_5_varscan$location) #2516
pat_10_liver_5_vs_fb <- intersect(pat_10_liver_5_varscan$location, pat_10_liver_5_freebayes$location) #2436

#intersect b/w mutect calls and two out of the 3 others
pat_10_liver_5_mutect_ug_fb <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_ug_fb) #6711
pat_10_liver_5_mutect_ug_vs <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_ug_vs) #528
pat_10_liver_5_mutect_vs_fb <- intersect(pat_10_liver_5_mutect$location, pat_10_liver_5_vs_fb) #495

#combine them
pat_10_liver_5_all <- c(pat_10_liver_5_mutect_ug_fb, pat_10_liver_5_mutect_ug_vs, pat_10_liver_5_mutect_vs_fb)

#subset mutect calls
pat_10_liver_5_mutect_short <- pat_10_liver_5_mutect[pat_10_liver_5_mutect$location %in% pat_10_liver_5_all, ] #6832

#trim to trusight gene calls
pat_10_liver_5_mutect_trusight <- trusight_from_mutect(pat_10_liver_5_mutect_short) #17

## plasma----
# import data
pat_10_plasma_ug <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #45442
pat_10_plasma_varscan <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #47419
pat_10_plasma_freebayes <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #356698

pat_10_plasma_mutect <- read.delim('/Volumes/Samsung_T5/pat_10/pat_10_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #575196

# dummy variable for matching
pat_10_plasma_ug$location <- paste0(pat_10_plasma_ug$CHROM, ':', pat_10_plasma_ug$POS)
pat_10_plasma_varscan$location <- paste0(pat_10_plasma_varscan$CHROM, ':', pat_10_plasma_varscan$POS)
pat_10_plasma_freebayes$location <- paste0(pat_10_plasma_freebayes$CHROM, ':', pat_10_plasma_freebayes$POS)

pat_10_plasma_mutect$location <- paste0(pat_10_plasma_mutect$CHROM, ':', pat_10_plasma_mutect$POS)

#pairwise intersections for the 3
pat_10_plasma_ug_fb <- intersect(pat_10_plasma_ug$location, pat_10_plasma_freebayes$location) #34782
pat_10_plasma_ug_vs <- intersect(pat_10_plasma_ug$location, pat_10_plasma_varscan$location) #5624
pat_10_plasma_vs_fb <- intersect(pat_10_plasma_varscan$location, pat_10_plasma_freebayes$location) #5250

#intersect b/w mutect calls and two out of the 3 others
pat_10_plasma_mutect_ug_fb <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_ug_fb) #25324
pat_10_plasma_mutect_ug_vs <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_ug_vs) #3249
pat_10_plasma_mutect_vs_fb <- intersect(pat_10_plasma_mutect$location, pat_10_plasma_vs_fb) #2915

#combine them
pat_10_plasma_all <- c(pat_10_plasma_mutect_ug_fb, pat_10_plasma_mutect_ug_vs, pat_10_plasma_mutect_vs_fb)

#subset mutect calls
pat_10_plasma_mutect_short <- pat_10_plasma_mutect[pat_10_plasma_mutect$location %in% pat_10_plasma_all, ] #25932

#trim to trusight gene calls
pat_10_plasma_mutect_trusight <- trusight_from_mutect(pat_10_plasma_mutect_short) #85

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 liver mets
pat_10_met_pool <- c(pat_10_liver_1_mutect_trusight$location, pat_10_liver_2a_mutect_trusight$location, pat_10_liver_5_mutect_trusight$location)
pat_10_met_pool <- unique(pat_10_met_pool) #23

pat_10_plasma_link <- pat_10_plasma_mutect_trusight[pat_10_plasma_mutect_trusight$location %in% pat_10_met_pool, ] #15 mutations in 9 genes









## PATIENT 9 ----
## lymph node----
# import data
pat_9_ln_ug <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ln_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4463
pat_9_ln_varscan <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ln_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #121
pat_9_ln_freebayes <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ln_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3975

pat_9_ln_mutect <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ln_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #2132

# dummy variable for matching
pat_9_ln_ug$location <- paste0(pat_9_ln_ug$CHROM, ':', pat_9_ln_ug$POS)
pat_9_ln_varscan$location <- paste0(pat_9_ln_varscan$CHROM, ':', pat_9_ln_varscan$POS)
pat_9_ln_freebayes$location <- paste0(pat_9_ln_freebayes$CHROM, ':', pat_9_ln_freebayes$POS)

pat_9_ln_mutect$location <- paste0(pat_9_ln_mutect$CHROM, ':', pat_9_ln_mutect$POS)

#pairwise intersections for the 3
pat_9_ln_ug_fb <- intersect(pat_9_ln_ug$location, pat_9_ln_freebayes$location) #1939
pat_9_ln_ug_vs <- intersect(pat_9_ln_ug$location, pat_9_ln_varscan$location) #66
pat_9_ln_vs_fb <- intersect(pat_9_ln_varscan$location, pat_9_ln_freebayes$location) #64

#intersect b/w mutect calls and two out of the 3 others
pat_9_ln_mutect_ug_fb <- intersect(pat_9_ln_mutect$location, pat_9_ln_ug_fb) #473
pat_9_ln_mutect_ug_vs <- intersect(pat_9_ln_mutect$location, pat_9_ln_ug_vs) #24
pat_9_ln_mutect_vs_fb <- intersect(pat_9_ln_mutect$location, pat_9_ln_vs_fb) #22

#combine them
pat_9_ln_all <- c(pat_9_ln_mutect_ug_fb, pat_9_ln_mutect_ug_vs, pat_9_ln_mutect_vs_fb)

#subset mutect calls
pat_9_ln_mutect_short <- pat_9_ln_mutect[pat_9_ln_mutect$location %in% pat_9_ln_all, ] #475

#trim to trusight gene calls
pat_9_ln_mutect_trusight <- trusight_from_mutect(pat_9_ln_mutect_short) #10

## oment----
# import data
pat_9_oment_ug <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_oment_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11359
pat_9_oment_varscan <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_oment_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #4839
pat_9_oment_freebayes <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_oment_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #88156

pat_9_oment_mutect <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_oment_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #114997

# dummy variable for matching
pat_9_oment_ug$location <- paste0(pat_9_oment_ug$CHROM, ':', pat_9_oment_ug$POS)
pat_9_oment_varscan$location <- paste0(pat_9_oment_varscan$CHROM, ':', pat_9_oment_varscan$POS)
pat_9_oment_freebayes$location <- paste0(pat_9_oment_freebayes$CHROM, ':', pat_9_oment_freebayes$POS)

pat_9_oment_mutect$location <- paste0(pat_9_oment_mutect$CHROM, ':', pat_9_oment_mutect$POS)

#pairwise intersections for the 3
pat_9_oment_ug_fb <- intersect(pat_9_oment_ug$location, pat_9_oment_freebayes$location) #8155
pat_9_oment_ug_vs <- intersect(pat_9_oment_ug$location, pat_9_oment_varscan$location) #1761
pat_9_oment_vs_fb <- intersect(pat_9_oment_varscan$location, pat_9_oment_freebayes$location) #1643

#intersect b/w mutect calls and two out of the 3 others
pat_9_oment_mutect_ug_fb <- intersect(pat_9_oment_mutect$location, pat_9_oment_ug_fb) #3432
pat_9_oment_mutect_ug_vs <- intersect(pat_9_oment_mutect$location, pat_9_oment_ug_vs) #575
pat_9_oment_mutect_vs_fb <- intersect(pat_9_oment_mutect$location, pat_9_oment_vs_fb) #507

#combine them
pat_9_oment_all <- c(pat_9_oment_mutect_ug_fb, pat_9_oment_mutect_ug_vs, pat_9_oment_mutect_vs_fb)

#subset mutect calls
pat_9_oment_mutect_short <- pat_9_oment_mutect[pat_9_oment_mutect$location %in% pat_9_oment_all, ] #3536

#trim to trusight gene calls
pat_9_oment_mutect_trusight <- trusight_from_mutect(pat_9_oment_mutect_short) #40

## ovary----
# import data
pat_9_ovary_ug <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ovary_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #8037
pat_9_ovary_varscan <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ovary_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #7513
pat_9_ovary_freebayes <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ovary_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #96219

pat_9_ovary_mutect <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_ovary_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #92251

# dummy variable for matching
pat_9_ovary_ug$location <- paste0(pat_9_ovary_ug$CHROM, ':', pat_9_ovary_ug$POS)
pat_9_ovary_varscan$location <- paste0(pat_9_ovary_varscan$CHROM, ':', pat_9_ovary_varscan$POS)
pat_9_ovary_freebayes$location <- paste0(pat_9_ovary_freebayes$CHROM, ':', pat_9_ovary_freebayes$POS)

pat_9_ovary_mutect$location <- paste0(pat_9_ovary_mutect$CHROM, ':', pat_9_ovary_mutect$POS)

#pairwise intersections for the 3
pat_9_ovary_ug_fb <- intersect(pat_9_ovary_ug$location, pat_9_ovary_freebayes$location) #6580
pat_9_ovary_ug_vs <- intersect(pat_9_ovary_ug$location, pat_9_ovary_varscan$location) #1391
pat_9_ovary_vs_fb <- intersect(pat_9_ovary_varscan$location, pat_9_ovary_freebayes$location) #1300

#intersect b/w mutect calls and two out of the 3 others
pat_9_ovary_mutect_ug_fb <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_ug_fb) #2981
pat_9_ovary_mutect_ug_vs <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_ug_vs) #519
pat_9_ovary_mutect_vs_fb <- intersect(pat_9_ovary_mutect$location, pat_9_ovary_vs_fb) #465

#combine them
pat_9_ovary_all <- c(pat_9_ovary_mutect_ug_fb, pat_9_ovary_mutect_ug_vs, pat_9_ovary_mutect_vs_fb)

#subset mutect calls
pat_9_ovary_mutect_short <- pat_9_ovary_mutect[pat_9_ovary_mutect$location %in% pat_9_ovary_all, ] #3055

#trim to trusight gene calls
pat_9_ovary_mutect_trusight <- trusight_from_mutect(pat_9_ovary_mutect_short) #20

## plasma----
# import data
pat_9_plasma_ug <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #59407
pat_9_plasma_varscan <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #20319
pat_9_plasma_freebayes <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #427734

pat_9_plasma_mutect <- read.delim('/Volumes/Samsung_T5/pat_9/pat_9_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #782680

# dummy variable for matching
pat_9_plasma_ug$location <- paste0(pat_9_plasma_ug$CHROM, ':', pat_9_plasma_ug$POS)
pat_9_plasma_varscan$location <- paste0(pat_9_plasma_varscan$CHROM, ':', pat_9_plasma_varscan$POS)
pat_9_plasma_freebayes$location <- paste0(pat_9_plasma_freebayes$CHROM, ':', pat_9_plasma_freebayes$POS)

pat_9_plasma_mutect$location <- paste0(pat_9_plasma_mutect$CHROM, ':', pat_9_plasma_mutect$POS)

#pairwise intersections for the 3
pat_9_plasma_ug_fb <- intersect(pat_9_plasma_ug$location, pat_9_plasma_freebayes$location) #52862
pat_9_plasma_ug_vs <- intersect(pat_9_plasma_ug$location, pat_9_plasma_varscan$location) #6420
pat_9_plasma_vs_fb <- intersect(pat_9_plasma_varscan$location, pat_9_plasma_freebayes$location) #6347

#intersect b/w mutect calls and two out of the 3 others
pat_9_plasma_mutect_ug_fb <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_ug_fb) #40548
pat_9_plasma_mutect_ug_vs <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_ug_vs) #2990
pat_9_plasma_mutect_vs_fb <- intersect(pat_9_plasma_mutect$location, pat_9_plasma_vs_fb) #2770

#combine them
pat_9_plasma_all <- c(pat_9_plasma_mutect_ug_fb, pat_9_plasma_mutect_ug_vs, pat_9_plasma_mutect_vs_fb)

#subset mutect calls
pat_9_plasma_mutect_short <- pat_9_plasma_mutect[pat_9_plasma_mutect$location %in% pat_9_plasma_all, ] #40970

#trim to trusight gene calls
pat_9_plasma_mutect_trusight <- trusight_from_mutect(pat_9_plasma_mutect_short) #59

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_9_met_pool <- c(pat_9_ln_mutect_trusight$location, pat_9_oment_mutect_trusight$location, pat_9_ovary_mutect_trusight$location)
pat_9_met_pool <- unique(pat_9_met_pool) #58

pat_9_plasma_link <- pat_9_plasma_mutect_trusight[pat_9_plasma_mutect_trusight$location %in% pat_9_met_pool, ] #32 mutations in 11 genes






## PATIENT 8 ----
## axillary----
# import data
pat_8_axillary_ug <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_axillary_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #7582
pat_8_axillary_varscan <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_axillary_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3213
pat_8_axillary_freebayes <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_axillary_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #54045

pat_8_axillary_mutect <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_axillary_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #68044

# dummy variable for matching
pat_8_axillary_ug$location <- paste0(pat_8_axillary_ug$CHROM, ':', pat_8_axillary_ug$POS)
pat_8_axillary_varscan$location <- paste0(pat_8_axillary_varscan$CHROM, ':', pat_8_axillary_varscan$POS)
pat_8_axillary_freebayes$location <- paste0(pat_8_axillary_freebayes$CHROM, ':', pat_8_axillary_freebayes$POS)

pat_8_axillary_mutect$location <- paste0(pat_8_axillary_mutect$CHROM, ':', pat_8_axillary_mutect$POS)

#pairwise intersections for the 3
pat_8_axillary_ug_fb <- intersect(pat_8_axillary_ug$location, pat_8_axillary_freebayes$location) #6239
pat_8_axillary_ug_vs <- intersect(pat_8_axillary_ug$location, pat_8_axillary_varscan$location) #1773
pat_8_axillary_vs_fb <- intersect(pat_8_axillary_varscan$location, pat_8_axillary_freebayes$location) #22

#intersect b/w mutect calls and two out of the 3 others
pat_8_axillary_mutect_ug_fb <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_ug_fb) #1937
pat_8_axillary_mutect_ug_vs <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_ug_vs) #271
pat_8_axillary_mutect_vs_fb <- intersect(pat_8_axillary_mutect$location, pat_8_axillary_vs_fb) #235

#combine them
pat_8_axillary_all <- c(pat_8_axillary_mutect_ug_fb, pat_8_axillary_mutect_ug_vs, pat_8_axillary_mutect_vs_fb)

#subset mutect calls
pat_8_axillary_mutect_short <- pat_8_axillary_mutect[pat_8_axillary_mutect$location %in% pat_8_axillary_all, ] #1987

#trim to trusight gene calls
pat_8_axillary_mutect_trusight <- trusight_from_mutect(pat_8_axillary_mutect_short) #12

## breast 1----
# import data
pat_8_breast_1_ug <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_1_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #11948
pat_8_breast_1_varscan <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_1_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3669
pat_8_breast_1_freebayes <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_1_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #74674

pat_8_breast_1_mutect <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_1_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #76214

# dummy variable for matching
pat_8_breast_1_ug$location <- paste0(pat_8_breast_1_ug$CHROM, ':', pat_8_breast_1_ug$POS)
pat_8_breast_1_varscan$location <- paste0(pat_8_breast_1_varscan$CHROM, ':', pat_8_breast_1_varscan$POS)
pat_8_breast_1_freebayes$location <- paste0(pat_8_breast_1_freebayes$CHROM, ':', pat_8_breast_1_freebayes$POS)

pat_8_breast_1_mutect$location <- paste0(pat_8_breast_1_mutect$CHROM, ':', pat_8_breast_1_mutect$POS)

#pairwise intersections for the 3
pat_8_breast_1_ug_fb <- intersect(pat_8_breast_1_ug$location, pat_8_breast_1_freebayes$location) #8351
pat_8_breast_1_ug_vs <- intersect(pat_8_breast_1_ug$location, pat_8_breast_1_varscan$location) #1759
pat_8_breast_1_vs_fb <- intersect(pat_8_breast_1_varscan$location, pat_8_breast_1_freebayes$location) #1774

#intersect b/w mutect calls and two out of the 3 others
pat_8_breast_1_mutect_ug_fb <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_ug_fb) #3080
pat_8_breast_1_mutect_ug_vs <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_ug_vs) #502
pat_8_breast_1_mutect_vs_fb <- intersect(pat_8_breast_1_mutect$location, pat_8_breast_1_vs_fb) #470

#combine them
pat_8_breast_1_all <- c(pat_8_breast_1_mutect_ug_fb, pat_8_breast_1_mutect_ug_vs, pat_8_breast_1_mutect_vs_fb)

#subset mutect calls
pat_8_breast_1_mutect_short <- pat_8_breast_1_mutect[pat_8_breast_1_mutect$location %in% pat_8_breast_1_all, ] #3136

#trim to trusight gene calls
pat_8_breast_1_mutect_trusight <- trusight_from_mutect(pat_8_breast_1_mutect_short) #19

## breast 2----
# import data
pat_8_breast_2_ug <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_2_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #62413
pat_8_breast_2_varscan <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_2_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3669
pat_8_breast_2_freebayes <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_2_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #70717

pat_8_breast_2_mutect <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_breast_2_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #61361

# dummy variable for matching
pat_8_breast_2_ug$location <- paste0(pat_8_breast_2_ug$CHROM, ':', pat_8_breast_2_ug$POS)
pat_8_breast_2_varscan$location <- paste0(pat_8_breast_2_varscan$CHROM, ':', pat_8_breast_2_varscan$POS)
pat_8_breast_2_freebayes$location <- paste0(pat_8_breast_2_freebayes$CHROM, ':', pat_8_breast_2_freebayes$POS)

pat_8_breast_2_mutect$location <- paste0(pat_8_breast_2_mutect$CHROM, ':', pat_8_breast_2_mutect$POS)

#pairwise intersections for the 3
pat_8_breast_2_ug_fb <- intersect(pat_8_breast_2_ug$location, pat_8_breast_2_freebayes$location) #9882
pat_8_breast_2_ug_vs <- intersect(pat_8_breast_2_ug$location, pat_8_breast_2_varscan$location) #1126
pat_8_breast_2_vs_fb <- intersect(pat_8_breast_2_varscan$location, pat_8_breast_2_freebayes$location) #1216

#intersect b/w mutect calls and two out of the 3 others
pat_8_breast_2_mutect_ug_fb <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_ug_fb) #2995
pat_8_breast_2_mutect_ug_vs <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_ug_vs) #850
pat_8_breast_2_mutect_vs_fb <- intersect(pat_8_breast_2_mutect$location, pat_8_breast_2_vs_fb) #795

#combine them
pat_8_breast_2_all <- c(pat_8_breast_2_mutect_ug_fb, pat_8_breast_2_mutect_ug_vs, pat_8_breast_2_mutect_vs_fb)

#subset mutect calls
pat_8_breast_2_mutect_short <- pat_8_breast_2_mutect[pat_8_breast_2_mutect$location %in% pat_8_breast_2_all, ] #3064

#trim to trusight gene calls
pat_8_breast_2_mutect_trusight <- trusight_from_mutect(pat_8_breast_2_mutect_short) #54

## plasma----
# import data
pat_8_plasma_ug <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_plasma_ann_ug_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #5277
pat_8_plasma_varscan <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_plasma_ann_varscan_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #1837
pat_8_plasma_freebayes <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_plasma_ann_freebayes_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #3975

pat_8_plasma_mutect <- read.delim('/Volumes/Samsung_T5/pat_8/pat_8_plasma_ann_mutect_vcf.txt', header = TRUE, stringsAsFactors = FALSE) #15390

# dummy variable for matching
pat_8_plasma_ug$location <- paste0(pat_8_plasma_ug$CHROM, ':', pat_8_plasma_ug$POS)
pat_8_plasma_varscan$location <- paste0(pat_8_plasma_varscan$CHROM, ':', pat_8_plasma_varscan$POS)
pat_8_plasma_freebayes$location <- paste0(pat_8_plasma_freebayes$CHROM, ':', pat_8_plasma_freebayes$POS)

pat_8_plasma_mutect$location <- paste0(pat_8_plasma_mutect$CHROM, ':', pat_8_plasma_mutect$POS)

#pairwise intersections for the 3
pat_8_plasma_ug_fb <- intersect(pat_8_plasma_ug$location, pat_8_plasma_freebayes$location) #1343
pat_8_plasma_ug_vs <- intersect(pat_8_plasma_ug$location, pat_8_plasma_varscan$location) #437
pat_8_plasma_vs_fb <- intersect(pat_8_plasma_varscan$location, pat_8_plasma_freebayes$location) #464

#intersect b/w mutect calls and two out of the 3 others
pat_8_plasma_mutect_ug_fb <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_ug_fb) #3386
pat_8_plasma_mutect_ug_vs <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_ug_vs) #274
pat_8_plasma_mutect_vs_fb <- intersect(pat_8_plasma_mutect$location, pat_8_plasma_vs_fb) #258

#combine them
pat_8_plasma_all <- c(pat_8_plasma_mutect_ug_fb, pat_8_plasma_mutect_ug_vs, pat_8_plasma_mutect_vs_fb)

#subset mutect calls
pat_8_plasma_mutect_short <- pat_8_plasma_mutect[pat_8_plasma_mutect$location %in% pat_8_plasma_all, ] #1367

#trim to trusight gene calls
pat_8_plasma_mutect_trusight <- trusight_from_mutect(pat_8_plasma_mutect_short) #13

## looking at how well plasma detects tumor mutations ----

#pool mutations from 3 mets
pat_8_met_pool <- c(pat_8_axillary_mutect_trusight$location, pat_8_breast_1_mutect_trusight$location, pat_8_breast_2_mutect_trusight$location)
pat_8_met_pool <- unique(pat_8_met_pool) #65

pat_8_plasma_link <- pat_8_plasma_mutect_trusight[pat_8_plasma_mutect_trusight$location %in% pat_8_met_pool, ] #10 mutations in 8 genes

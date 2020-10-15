# mutation analysis
# bring in data first
# there are 4 patients, 010, 008, 002, and 009 (which has 2 sets of data)

# for each patient, comparisons done:
#   any tumor (primary or met) vs buffy coat
#   plasma vs buffy coat

# do mets/tumors have similar mutations?
# are any also similar to plasma?

##libraries----
library(readxl)
library(ggsci)
library(UpSetR)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(gplots)
library(plotrix)
library(stringr)

##functions for processing mutation files----
mutation_process_1 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(20:23)]
  mutation_file <- mutation_file[mutation_file$Function != 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  mutation_file <- mutation_file[mutation_file$chrom != 'chrY', ]
  return(mutation_file)
}

mutation_process_2 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(19:22)]
  mutation_file <- mutation_file[mutation_file$Func != 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  mutation_file <- mutation_file[mutation_file$chrom != 'chrY', ]
  return(mutation_file)
}

mutation_process_3 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(19:22)]
  mutation_file <- mutation_file[mutation_file$func != 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  mutation_file <- mutation_file[mutation_file$chrom != 'chrY', ]
  return(mutation_file)
}
## patient 10 data----
# this patient has 3 liver mets and plasma
pat_10_liv_met_2a <- mutation_process_1('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP9_snv_annot.xlsx', 1)
# 246 left

pat_10_liv_met_1 <- mutation_process_1('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP10_snv_annot.xlsx', 1)
# 233 left

pat_10_liv_met_5 <- mutation_process_2('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP12_snv_annot.xlsx', 1)
# 254 left

pat_10_plasma <- mutation_process_2('N:/D13100 BGI analysis/Pt 010/CP2-VS-UP1_snv_annot (low AF).xlsx', 1)
# 352 left 96027 using low AF

#what do they have in common
all_pat_10_genes <- Reduce(intersect, list(pat_10_liv_met_1$`Transcript change`, pat_10_liv_met_2a$`Transcript change`, pat_10_liv_met_5$Transcript))
#3

#are these also in plasma?
plasma_intersect_pat_10 <- intersect(all_pat_10_genes, pat_10_plasma$Transcript)
#all 3 are there


# subset to these genes and take a look
pat_10_liv_met_1_short <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% plasma_intersect_pat_10, ]

pat_10_liv_met_2a_short <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% plasma_intersect_pat_10, ]

pat_10_liv_met_5_short <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% plasma_intersect_pat_10, ]

pat_10_plasma_short <- pat_10_plasma[pat_10_plasma$Transcript %in% plasma_intersect_pat_10, ]

# also compare plasma to pool of tumors
pat_10_pooled <- unique(c(pat_10_liv_met_1$`Transcript change`, pat_10_liv_met_2a$`Transcript change`, pat_10_liv_met_5$Transcript)) #712 unique mutations

pat_10_pooled_plasma <- intersect(pat_10_pooled, pat_10_plasma$Transcript) #27 or 81 with low AF

# subset to pooled plasma and take a look
pat_10_liv_met_1_pooled <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_pooled_plasma, ]
#12 or 32
pat_10_liv_met_2a_pooled <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_pooled_plasma, ]
#8 or 29
pat_10_liv_met_5_pooled <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_pooled_plasma, ]
#15 or 33
pat_10_plasma_pooled <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_pooled_plasma, ]
#27 or 81

# subset these to tumor AF for new figure
#pat_10_liv_met_1_pooled <- pat_10_liv_met_1_pooled[pat_10_liv_met_1_pooled$Function == 'missense_variant', ]
pat_10_liv_met_1_pooled <- pat_10_liv_met_1_pooled[, c('Gene', 'Transcript change', 'tumor_var_freq')]
pat_10_liv_met_1_pooled$`Transcript change` <- gsub('^.*:', '', pat_10_liv_met_1_pooled$`Transcript change`)
pat_10_liv_met_1_pooled$`Transcript change` <- gsub('\\/.*$', '', pat_10_liv_met_1_pooled$`Transcript change`)
pat_10_liv_met_1_pooled$`Transcript change` <- paste0(pat_10_liv_met_1_pooled$Gene, '\n', pat_10_liv_met_1_pooled$`Transcript change`)
# pat_10_liv_met_1_pooled$`Transcript change` <- str_sub(pat_10_liv_met_1_pooled$`Transcript change`, start = -3)
# pat_10_liv_met_1_pooled$`Transcript change` <- gsub('>', ':', pat_10_liv_met_1_pooled$`Transcript change`)
# pat_10_liv_met_1_pooled$`Transcript change` <- paste0(pat_10_liv_met_1_pooled$chrom, ':', pat_10_liv_met_1_pooled$position, ':', pat_10_liv_met_1_pooled$`Transcript change`)
pat_10_liv_met_1_pooled <- pat_10_liv_met_1_pooled[, c('Transcript change', 'tumor_var_freq')]
colnames(pat_10_liv_met_1_pooled) <- c('Transcript', 'pat_10_liv_met_1')

#pat_10_liv_met_2a_pooled <- pat_10_liv_met_2a_pooled[pat_10_liv_met_2a_pooled$Function == 'missense_variant', ]
pat_10_liv_met_2a_pooled <- pat_10_liv_met_2a_pooled[, c('Gene', 'Transcript change', 'tumor_var_freq')]
pat_10_liv_met_2a_pooled$`Transcript change` <- gsub('^.*:', '', pat_10_liv_met_2a_pooled$`Transcript change`)
pat_10_liv_met_2a_pooled$`Transcript change` <- gsub('\\/.*$', '', pat_10_liv_met_2a_pooled$`Transcript change`)
pat_10_liv_met_2a_pooled$`Transcript change` <- paste0(pat_10_liv_met_2a_pooled$Gene, '\n', pat_10_liv_met_2a_pooled$`Transcript change`)
# pat_10_liv_met_2a_pooled$`Transcript change` <- str_sub(pat_10_liv_met_2a_pooled$`Transcript change`, start = -3)
# pat_10_liv_met_2a_pooled$`Transcript change` <- gsub('>', ':', pat_10_liv_met_2a_pooled$`Transcript change`)
# pat_10_liv_met_2a_pooled$`Transcript change` <- paste0(pat_10_liv_met_2a_pooled$chrom, ':', pat_10_liv_met_2a_pooled$position, ':', pat_10_liv_met_2a_pooled$`Transcript change`)
pat_10_liv_met_2a_pooled <- pat_10_liv_met_2a_pooled[, c('Transcript change', 'tumor_var_freq')]
colnames(pat_10_liv_met_2a_pooled) <- c('Transcript', 'pat_10_liv_met_2a')

#pat_10_liv_met_5_pooled <- pat_10_liv_met_5_pooled[pat_10_liv_met_5_pooled$Func == 'missense_variant', ]
pat_10_liv_met_5_pooled <- pat_10_liv_met_5_pooled[, c('Gene', 'Transcript', 'tumor_var_freq')]
pat_10_liv_met_5_pooled$Transcript <- gsub('^.*:', '', pat_10_liv_met_5_pooled$Transcript)
pat_10_liv_met_5_pooled$Transcript <- gsub('\\/.*$', '', pat_10_liv_met_5_pooled$Transcript)
pat_10_liv_met_5_pooled$Transcript <- paste0(pat_10_liv_met_5_pooled$Gene, '\n', pat_10_liv_met_5_pooled$Transcript)
# pat_10_liv_met_5_pooled$Transcript <- str_sub(pat_10_liv_met_5_pooled$Transcript, start = -3)
# pat_10_liv_met_5_pooled$Transcript <- gsub('>', ':', pat_10_liv_met_5_pooled$Transcript)
# pat_10_liv_met_5_pooled$Transcript <- paste0(pat_10_liv_met_5_pooled$chrom, ':', pat_10_liv_met_5_pooled$position, ':', pat_10_liv_met_5_pooled$Transcript)
pat_10_liv_met_5_pooled <- pat_10_liv_met_5_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_10_liv_met_5_pooled) <- c('Transcript', 'pat_10_liv_met_5')

#pat_10_plasma_pooled <- pat_10_plasma_pooled[pat_10_plasma_pooled$Func == 'missense_variant', ]
pat_10_plasma_pooled <- pat_10_plasma_pooled[, c('Gene', 'Transcript', 'tumor_var_freq')]
pat_10_plasma_pooled$Transcript <- gsub('^.*:', '', pat_10_plasma_pooled$Transcript)
pat_10_plasma_pooled$Transcript <- gsub('\\/.*$', '', pat_10_plasma_pooled$Transcript)
pat_10_plasma_pooled$Transcript <- paste0(pat_10_plasma_pooled$Gene, '\n', pat_10_plasma_pooled$Transcript)
# pat_10_plasma_pooled$Transcript <- str_sub(pat_10_plasma_pooled$Transcript, start = -3)
# pat_10_plasma_pooled$Transcript <- gsub('>', ':', pat_10_plasma_pooled$Transcript)
# pat_10_plasma_pooled$Transcript <- paste0(pat_10_plasma_pooled$chrom, ':', pat_10_plasma_pooled$position, ':', pat_10_plasma_pooled$Transcript)
pat_10_plasma_pooled <- pat_10_plasma_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_10_plasma_pooled) <- c('Transcript', 'pat_10_plasma')

# put them together
pat_10_muts_pooled <- merge(pat_10_liv_met_1_pooled, pat_10_liv_met_2a_pooled, by = 'Transcript', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_liv_met_5_pooled, by = 'Transcript', all = TRUE)
pat_10_muts_pooled <- merge(pat_10_muts_pooled, pat_10_plasma_pooled, by = 'Transcript', all = TRUE)
pat_10_muts_pooled <- pat_10_muts_pooled[order(-pat_10_muts_pooled$pat_10_plasma), ]

rownames(pat_10_muts_pooled) <- pat_10_muts_pooled$Transcript
row_muts <- rownames(pat_10_muts_pooled)
pat_10_muts_pooled <- pat_10_muts_pooled[, -1]
#pat_10_muts_pooled <- sapply(pat_10_muts_pooled, function(x) ifelse (is.na(x), 0, x))
pat_10_muts_pooled <- as.data.frame(pat_10_muts_pooled)
rownames(pat_10_muts_pooled) <- row_muts
#plot figure to replace in manuscript

# 2d heatmap
# subset out plasma col
pat_10_3 <- c(pat_10_muts_pooled$pat_10_liv_met_1, pat_10_muts_pooled$pat_10_liv_met_2a, pat_10_muts_pooled$pat_10_liv_met_5)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_10_muts_pooled)[1] * dim(pat_10_muts_pooled)[2])
# plasma reds
cell_cols[244:324] <- color.scale(pat_10_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:243] <- color.scale(pat_10_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 81, byrow = FALSE)
pat_10_pooled_t <- data.frame(t(pat_10_muts_pooled))
pat_10_pooled_t <- pat_10_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_10_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
legval<-seq(min(pat_10_muts_pooled[, 4], na.rm = TRUE),max(pat_10_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(0,-0.9,30,-0.5,round(c(min(pat_10_muts_pooled[, 4], na.rm = TRUE), max(pat_10_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.2, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE),max(pat_10_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(35,-0.9,65,-0.5,round(c(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE), max(pat_10_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=37, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 32.5, cex = 1.1, font = 2)

# add NA legend
color.legend(68, -0.9, 70, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=3.0, at=72.9, cex = 1.1, font = 2)
legend(x=68.15,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.9, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver\nMet 1', 'Liver\nMet 2', 'Liver\nMet 5')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_1))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_2a)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_2a))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_5)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_5))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))







## intergenic?----
intergenic_process_1 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(20:23)]
  mutation_file <- mutation_file[mutation_file$Function == 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  return(mutation_file)
}

intergenic_process_2 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(19:22)]
  mutation_file <- mutation_file[mutation_file$Func == 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  return(mutation_file)
}

intergenic_process_3 <- function(mutations_file, sheet) {
  mutation_file <- read_excel(mutations_file, sheet = sheet)
  mutation_file <- mutation_file[, -c(19:22)]
  mutation_file <- mutation_file[mutation_file$func == 'intergenic_region', ]
  mutation_file <- mutation_file[mutation_file$normal_var_freq < 0.01, ]
  var_idx <- mutation_file$tumor_reads2 == 2 & mutation_file$tumor_var_freq >= 0.165
  mutation_file <- mutation_file[!var_idx, ]
  return(mutation_file)
}

pat_10_liv_met_2a_intergenic <- intergenic_process_1('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP9_snv_annot.xlsx', 1)
# 1136 left

pat_10_liv_met_1_intergenic <- intergenic_process_1('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP10_snv_annot.xlsx', 1)
# 1365 left

pat_10_liv_met_5_intergenic <- intergenic_process_2('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP12_snv_annot.xlsx', 1)
# 1231 left

pat_10_plasma_intergenic <- intergenic_process_2('N:/D13100 BGI analysis/Pt 010/CP2-VS-UP1_snv_annot.xlsx', 1)
# 1446 left

# what do they have in common?
all_pat_10_genes_intergenic <- Reduce(intersect, list(pat_10_liv_met_1_intergenic$`Transcript change`, pat_10_liv_met_2a_intergenic$`Transcript change`, 
                                                      pat_10_liv_met_5_intergenic$Transcript)) #21

#are these also in plasma?
plasma_intersect_pat_10_intergenic <- intersect(all_pat_10_genes, pat_10_plasma_intergenic$Transcript)
#NONE

# subset to these genes and take a look
pat_10_liv_met_1_short_intergenic <- pat_10_liv_met_1_intergenic[pat_10_liv_met_1_intergenic$`Transcript change` %in% all_pat_10_genes_intergenic, ]
#21

pat_10_liv_met_2a_short_intergenic <- pat_10_liv_met_2a_intergenic[pat_10_liv_met_2a_intergenic$`Transcript change` %in% all_pat_10_genes_intergenic, ]
#21

pat_10_liv_met_5_short_intergenic <- pat_10_liv_met_5_intergenic[pat_10_liv_met_5_intergenic$Transcript %in% all_pat_10_genes_intergenic, ]
#21

#pat_10_plasma_short <- pat_10_plasma[pat_10_plasma$Transcript %in% plasma_intersect_pat_10, ]

##function histograms----
Intergenic <- rep(c('Intergenic', 'Non-Intergenic'), 4)
Sample <- c(rep('Liver_Met_1', 2), rep('Liver_Met_2a', 2), rep('Liver_Met_5', 2), rep('Plasma', 2))
value <- c(1136, 234, 1365, 247, 1231, 255, 1446, 353)
dfa <- data.frame(Intergenic, Sample, value)
p <- ggplot(dfa) + scale_y_continuous(limits = c(0, 1799), breaks = NULL)
p1 <- p + geom_bar(aes(Sample, value), stat = "identity", fill = 'purple')
p2 <- p1 + geom_bar(aes(Sample, value, fill = Intergenic),
                    stat = "identity", position = "dodge")
p2 +scale_fill_brewer(palette = 'Purples') + xlab('Samples') + ggtitle('Patient 10') + ylab('Count') + labs(fill = 'Functions') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white"))
# do this for non-intergenic separately
pat_10_liv_met_1_df <- read_excel('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP10_snv_annot.xlsx')
pat_10_liv_met_1_df <- pat_10_liv_met_1_df[pat_10_liv_met_1_df$normal_var_freq < 0.01, ]
var_idx <- pat_10_liv_met_1_df$tumor_reads2 == 2 & pat_10_liv_met_1_df$tumor_var_freq >= 0.165
pat_10_liv_met_1_df <- pat_10_liv_met_1_df[!var_idx, ]
pat_10_classes <- table(pat_10_liv_met_1_df$Function)
pat_10_functions <- names(pat_10_classes)
pat_10_values <- as.vector(pat_10_classes)
pat_10_samples <- rep('liv_met_1', 10)
dfa <- data.frame(pat_10_functions, pat_10_samples, pat_10_values)

pat_10_liv_met_2_df <- read_excel('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP9_snv_annot.xlsx')
pat_10_liv_met_2_df <- pat_10_liv_met_2_df[pat_10_liv_met_2_df$normal_var_freq < 0.01, ]
var_idx <- pat_10_liv_met_2_df$tumor_reads2 == 2 & pat_10_liv_met_2_df$tumor_var_freq >= 0.165
pat_10_liv_met_2_df <- pat_10_liv_met_2_df[!var_idx, ]
pat_10_classes <- table(pat_10_liv_met_2_df$Function)
pat_10_functions <- names(pat_10_classes)
pat_10_values <- as.vector(pat_10_classes)
pat_10_samples <- rep('liv_met_2a', 11)
dfa2 <- data.frame(pat_10_functions, pat_10_samples, pat_10_values)

pat_10_liv_met_5_df <- read_excel('N:/D13100 BGI analysis/Pt 010/CP2-VS-CP12_snv_annot.xlsx')
pat_10_liv_met_5_df <- pat_10_liv_met_5_df[pat_10_liv_met_5_df$normal_var_freq < 0.01, ]
var_idx <- pat_10_liv_met_5_df$tumor_reads2 == 2 & pat_10_liv_met_5_df$tumor_var_freq >= 0.165
pat_10_liv_met_5_df <- pat_10_liv_met_5_df[!var_idx, ]
pat_10_classes <- table(pat_10_liv_met_5_df$Func)
pat_10_functions <- names(pat_10_classes)
pat_10_values <- as.vector(pat_10_classes)
pat_10_samples <- rep('liv_met_5', 11)
dfa3 <- data.frame(pat_10_functions, pat_10_samples, pat_10_values)

pat_10_plasma_df <- read_excel('N:/D13100 BGI analysis/Pt 010/CP2-VS-UP1_snv_annot.xlsx')
pat_10_plasma_df <- pat_10_plasma_df[pat_10_plasma_df$normal_var_freq < 0.01, ]
var_idx <- pat_10_plasma_df$tumor_reads2 == 2 & pat_10_plasma_df$tumor_var_freq >= 0.165
pat_10_plasma_df <- pat_10_plasma_df[!var_idx, ]
pat_10_classes <- table(pat_10_plasma_df$Func)
pat_10_functions <- names(pat_10_classes)
pat_10_values <- as.vector(pat_10_classes)
pat_10_samples <- rep('plasma', 12)
dfa4 <- data.frame(pat_10_functions, pat_10_samples, pat_10_values)

dfa <- rbind(dfa, dfa2)
dfa <- rbind(dfa, dfa3)
dfa <- rbind(dfa, dfa4)
dfa <- dfa[dfa$pat_10_functions != 'intergenic_region', ]

ylim <- max(cast(dfa, pat_10_samples ~ ., sum)[, 2],
            cast(dfa, pat_10_functions ~ ., sum)[, 2])
my_new_palette <- brewer.pal(11, 'Spectral')
my_new_palette <- colorRampPalette(my_new_palette)(13)

p <- ggplot(dfa) + scale_y_continuous(limits = c(0, 375), breaks = NULL)
p1 <- p + geom_bar(aes(pat_10_samples, pat_10_values), stat = "identity", fill = 'dodgerblue')
p2 <- p1 + geom_bar(aes(pat_10_samples, pat_10_values, fill = pat_10_functions),
                    stat = "identity", position = "dodge")
p2 +scale_fill_manual(values = my_new_palette) + xlab('Samples') + ggtitle('Patient 10') + ylab('Count') + labs(fill = 'Functions') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "white"))








## patient 2----
# this patient has 2 breast tumor samples (should be the same) and 2 liver mets plus plasma
# check these first for similarity
pat_2_breast_1 <- mutation_process_1('N:/D13100 BGI analysis/Pt 002/UP13-VS-UP9_snv_annot.xlsx', 1)
#236 left

pat_2_breast_2 <- mutation_process_1('N:/D13100 BGI analysis/Pt 002/UP13-VS-UP10_snv_annot.xlsx', 1)
#209 left

# mutations in common?
length(intersect(pat_2_breast_1$`Transcript change`, pat_2_breast_2$`Transcript change`))
# only 21 in common, these are clearly different!!!

# work up the 2 liver mets and look for similarities
pat_2_liv_met_1 <- mutation_process_1('N:/D13100 BGI analysis/Pt 002/UP13-VS-UP12_snv_annot.xlsx', 1)
# 267 left

pat_2_liv_met_2 <- mutation_process_1('N:/D13100 BGI analysis/Pt 002/UP13-VS-UP11_snv_annot.xlsx', 1)
# 396 left

#liver met common mutations
pat_2_liv_common <- intersect(pat_2_liv_met_1$`Transcript change`, pat_2_liv_met_2$`Transcript change`)
#13 in common

#how similar are each breast sample to the liver mets?
breast_1_common <- intersect(pat_2_liv_common, pat_2_breast_1$`Transcript change`)
#6 of the 13
breast_2_common <- intersect(pat_2_liv_common, pat_2_breast_2$`Transcript change`)
#6 of the 13

#are these the same 6??
intersect(breast_1_common, breast_2_common)
#5 of 6 in common
# use these?

#subset both based on common liver mets
pat_2_breast_1_short <- pat_2_breast_1[pat_2_breast_1$`Transcript change` %in% pat_2_liv_common, ]
#6
pat_2_breast_2_short <- pat_2_breast_2[pat_2_breast_2$`Transcript change` %in% pat_2_liv_common, ]
#6

#which has more possibly confident calls?
sum(pat_2_breast_1_short$tumor_reads2 == 2) #0
sum(pat_2_breast_2_short$tumor_reads2 == 2) #0
# use breast 2 and breast 2 common b/c the one not in common is a little more significant

pat_2_plasma <- mutation_process_2('N:/D13100 BGI analysis/Pt 002/UP13-VS-UP8_snv_annot (low AF).xlsx', 1)
#668 left or 848

plasma_intersect_pat_2 <- intersect(pat_2_plasma$Transcript, breast_2_common)
#NONE in common

#subset and look
pat_2_breast_2_shorter <- pat_2_breast_2_short[pat_2_breast_2_short$`Transcript change` %in% breast_2_common, ]
#6
pat_2_liv_met_1_short <- pat_2_liv_met_1[pat_2_liv_met_1$`Transcript change` %in% breast_2_common, ]
#6
pat_2_liv_met_2_short <- pat_2_liv_met_2[pat_2_liv_met_2$`Transcript change` %in% breast_2_common, ]
#6

#none of these are in plasma
#pat_2_plasma_short <- pat_2_plasma[pat_2_plasma$`Transcript change` %in% breast_2_common, ]

# look at pooled vs plasma
pat_2_pooled <- unique(c(pat_2_breast_2$`Transcript change`, pat_2_liv_met_1$`Transcript change`, pat_2_liv_met_2$`Transcript change`)) #841 unique mutations

pat_2_pooled_plasma <- intersect(pat_2_pooled, pat_2_plasma$Transcript) #1

# subset to pooled plasma and take a look
pat_2_breast_2_pooled <- pat_2_breast_2[pat_2_breast_2$`Transcript change` %in% pat_2_pooled_plasma, ]
#0
pat_2_liv_met_1_pooled <- pat_2_liv_met_1[pat_2_liv_met_1$`Transcript change` %in% pat_2_pooled_plasma, ]
#1
pat_2_liv_met_2_pooled <- pat_2_liv_met_2[pat_2_liv_met_2$`Transcript change` %in% pat_2_pooled_plasma, ]
#0
pat_2_plasma_pooled <- pat_2_plasma[pat_2_plasma$Transcript %in% pat_2_pooled_plasma, ]
#1

# subset these to tumor AF for new figure
pat_2_liv_met_1_pooled <- pat_2_liv_met_1_pooled[, c('position', 'chrom', 'Transcript change', 'tumor_var_freq')]
pat_2_liv_met_1_pooled$`Transcript change` <- str_sub(pat_2_liv_met_1_pooled$`Transcript change`, start = -3)
pat_2_liv_met_1_pooled$`Transcript change` <- gsub('>', ':', pat_2_liv_met_1_pooled$`Transcript change`)
pat_2_liv_met_1_pooled$`Transcript change` <- paste0(pat_2_liv_met_1_pooled$chrom, ':', pat_2_liv_met_1_pooled$position, ':', pat_2_liv_met_1_pooled$`Transcript change`)
pat_2_liv_met_1_pooled <- pat_2_liv_met_1_pooled[, c('Transcript change', 'tumor_var_freq')]
colnames(pat_2_liv_met_1_pooled) <- c('Transcript', 'pat_2_liv_met_1')

pat_2_plasma_pooled <- pat_2_plasma_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
pat_2_plasma_pooled$Transcript <- str_sub(pat_2_plasma_pooled$Transcript, start = -3)
pat_2_plasma_pooled$Transcript <- gsub('>', ':', pat_2_plasma_pooled$Transcript)
pat_2_plasma_pooled$Transcript <- paste0(pat_2_plasma_pooled$chrom, ':', pat_2_plasma_pooled$position, ':', pat_2_plasma_pooled$Transcript)
pat_2_plasma_pooled <- pat_2_plasma_pooled[, c('Transcript change', 'tumor_var_freq')]
colnames(pat_2_plasma_pooled) <- c('Transcript', 'pat_2_plasma')

# put them together
pat_2_muts_pooled <- merge(pat_2_liv_met_1_pooled, pat_2_plasma_pooled, by = 'Transcript', all = TRUE)


rownames(pat_2_muts_pooled) <- pat_2_muts_pooled$Transcript
row_muts <- rownames(pat_2_muts_pooled)
pat_2_muts_pooled <- pat_2_muts_pooled[, -1]
#pat_10_muts_pooled <- sapply(pat_10_muts_pooled, function(x) ifelse (is.na(x), 0, x))
pat_2_muts_pooled <- as.data.frame(pat_2_muts_pooled)
rownames(pat_2_muts_pooled) <- row_muts
#plot figure to replace in manuscript

# 2d heatmap
# subset out plasma col
pat_2_3 <- pat_2_muts_pooled$pat_2_liv_met_1

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_2_muts_pooled)[1] * dim(pat_2_muts_pooled)[2])
# plasma reds
cell_cols[2] <- color.scale(pat_2_muts_pooled[, 2], extremes = c('white', 'red'), na.color = '#00ff00')
# tumor blues
cell_cols[1] <- color.scale(pat_2_3, extremes = c('white', 'blue'), na.color = '#00ff00')
cell_cols <- matrix(cell_cols, nrow = 1, byrow = FALSE)
pat_2_pooled_t <- data.frame(t(pat_2_muts_pooled))
pat_2_pooled_t <- pat_2_pooled_t[c(2,1), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(2,1), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))

color2D.matplot(data.frame(pat_2_pooled_t), cellcolors=cell_cols, xlab = '', ylab = '', border=NA, axes = FALSE)

# add plasma legend
legval<-seq(min(pat_10_muts_pooled[, 4], na.rm = TRUE),max(pat_10_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('white', 'red'))
color.legend(0,-0.9,10,-0.5,round(c(min(pat_10_muts_pooled[, 4], na.rm = TRUE), max(pat_10_muts_pooled[, 4], na.rm = TRUE)),1),rect.col=legcol)
mtext('Plasma Mutant Allele Frequency', side=1, line=1.5, at=12.8)

# add tumor legend
legval<-seq(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE),max(pat_10_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('white', 'blue'))
color.legend(0,-1.7,10,-1.3,round(c(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE), max(pat_10_muts_pooled[, 1:3], na.rm = TRUE)),1),rect.col=legcol)
mtext('Tumor Mutant Allele Frequency', side=1, line=3.7, at=12.7)

# add NA legend
color.legend(15.5, -1.7, 15.9, -1.3, legend = '', rect.col = '#00ff00')
mtext('Mutation Not Present', side=1, line=3.7, at=17.7)
legend(x=15.35,y=-1.15,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
mut_col_end <- str_sub(mut_col_labels, -3)
mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.6, las = 2)

mut_row_labels <- c('Plasma', 'Liver\nMet 1', 'Liver\nMet 2', 'Liver\nMet 5')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1, las = 1)

#add points for NA values
points(x = c(0.5, 0.5, 1.5, 1.5, 2.5, 2.5, 4.5, 4.5, 5.5, 5.5, 6.5, 6.5, 7.5, 7.5, 8.5, 8.5, 9.5, 9.5, 10.5, 10.5), 
       y = c(2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 1.5, 0.5, 2.5, 1.5, 2.5, 1.5, 1.5, 0.5, 1.5, 0.5,  1.5,  0.5), 
       pch = 16)
points(x = c(11.5, 12.5, 12.5, 14.5, 14.5, 15.5, 16.5, 16.5, 17.5, 17.5, 18.5, 18.5, 19.5, 19.5, 20.5, 20.5), 
       y = c( 0.5,  2.5,  0.5,  1.5,  0.5,  0.5,  2.5,  1.5,  2.5,  0.5,  2.5,  0.5,  2.5,  1.5,  2.5,  1.5), 
       pch = 16)
points(x = c(21.5, 21.5, 22.5, 22.5, 23.5, 23.5, 24.5, 24.5, 25.5, 25.5), 
       y = c( 2.5,  1.5,  2.5,  1.5,  2.5,  1.5,  1.5,  0.5,  1.5,  0.5), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))




## patient 8----
# this patient has 2 breast tumors, an axillary met and plasma
pat_8_breast_1 <- mutation_process_3('N:/D13100 BGI analysis/Pt 008/UP4-VS-UP6_snv_annot.xlsx', 1)
# 444 left

pat_8_breast_2 <- mutation_process_3('N:/D13100 BGI analysis/Pt 008/UP4-VS-UP7_snv_annot.xlsx', 1)
# 1602 left

pat_8_axillary <- mutation_process_3('N:/D13100 BGI analysis/Pt 008/UP4-VS-UP5_snv_annot.xlsx', 1)
# 371 left

pat_8_plasma <- mutation_process_2('N:/D13100 BGI analysis/Pt 008/UP4-VS-UP3_snv_annot (low AF).xlsx', 1)
# 281 left or 312

#what do they have in common?
all_pat_8_genes <- Reduce(intersect, list(pat_8_breast_1$transcript, pat_8_breast_2$transcript, pat_8_axillary$transcript))
# 3 in common

#also in plasma?
intersect(pat_8_plasma$Transcript, all_pat_8_genes)
#0 NONE in plasma

# subset to common
pat_8_breast_1_short <- pat_8_breast_1[pat_8_breast_1$transcript %in% all_pat_8_genes, ]
#3

pat_8_breast_2_short <- pat_8_breast_2[pat_8_breast_2$transcript %in% all_pat_8_genes, ]
#3

pat_8_axillary_short <- pat_8_axillary[pat_8_axillary$transcript %in% all_pat_8_genes, ]
#3

# look at pooled vs plasma
pat_8_pooled <- unique(c(pat_8_breast_1$transcript, pat_8_breast_2$transcript, pat_8_axillary$transcript)) #2378 unique mutations

pat_8_pooled_plasma <- intersect(pat_8_pooled, pat_8_plasma$Transcript) #NONE


## patient 9----
# this patient has lymph, omental, and ovary mets plus plasma
pat_9_lymph_met <- mutation_process_3('N:/D13100 BGI analysis/Pt 009/CP1-VS-CP13_snv_annot.xlsx', 1)
# 565 left

pat_9_oment_met <- mutation_process_3('N:/D13100 BGI analysis/Pt 009/CP1-VS-CP14_snv_annot.xlsx', 1)
# 489 left

pat_9_ovary_met <- mutation_process_3('N:/D13100 BGI analysis/Pt 009/CP1-VS-CP15_snv_annot.xlsx', 1)
# 516 left

pat_9_plasma <- mutation_process_2('N:/D13100 BGI analysis/Pt 009/CP1-VS-CP16_snv_annot (low AF).xlsx', 1)
# 322 left or 41124

# what do they have in common
all_pat_9_genes <- Reduce(intersect, list(pat_9_lymph_met$transcript, pat_9_oment_met$transcript, pat_9_ovary_met$transcript))
# only 2

# in plasma?
plasma_intersect_pat_9 <- intersect(pat_9_plasma$Transcript, all_pat_9_genes)
# both in plasma

# subset and take a look
pat_9_lymph_met_short <- pat_9_lymph_met[pat_9_lymph_met$transcript %in% plasma_intersect_pat_9, ]
#2

pat_9_oment_met_short <- pat_9_oment_met[pat_9_oment_met$transcript %in% plasma_intersect_pat_9, ]
#2

pat_9_ovary_met_short <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% plasma_intersect_pat_9, ]
#2

pat_9_plasma_short <- pat_9_plasma[pat_9_plasma$Transcript %in% plasma_intersect_pat_9, ]
#2

# pooled vs plasma
pat_9_pooled <- unique(c(pat_9_lymph_met$transcript, pat_9_oment_met$transcript, pat_9_ovary_met$transcript)) #1551

pat_9_pooled_plasma <- intersect(pat_9_pooled, pat_9_plasma$Transcript) #13

# subset to pooled plasma and take a look
pat_9_lymph_met_pooled <- pat_9_lymph_met[pat_9_lymph_met$transcript %in% pat_9_pooled_plasma, ]
#2
pat_9_oment_met_pooled <- pat_9_oment_met[pat_9_oment_met$transcript %in% pat_9_pooled_plasma, ]
#8
pat_9_ovary_met_pooled <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% pat_9_pooled_plasma, ]
#10
pat_9_plasma_pooled <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_pooled_plasma, ]
#13

# subset these to tumor AF for new figure
pat_9_lymph_met_pooled <- pat_9_lymph_met_pooled[, c('position', 'chrom', 'transcript', 'tumor_var_freq')]
pat_9_lymph_met_pooled$transcript <- str_sub(pat_9_lymph_met_pooled$transcript, start = -3)
pat_9_lymph_met_pooled$transcript <- gsub('>', ':', pat_9_lymph_met_pooled$transcript)
pat_9_lymph_met_pooled$transcript <- paste0(pat_9_lymph_met_pooled$chrom, ':', pat_9_lymph_met_pooled$position, ':', pat_9_lymph_met_pooled$transcript)
pat_9_lymph_met_pooled <- pat_9_lymph_met_pooled[, c('transcript', 'tumor_var_freq')]
colnames(pat_9_lymph_met_pooled) <- c('Transcript', 'pat_9_lymph_met')

pat_9_oment_met_pooled <- pat_9_oment_met_pooled[, c('position', 'chrom', 'transcript', 'tumor_var_freq')]
pat_9_oment_met_pooled$transcript <- str_sub(pat_9_oment_met_pooled$transcript, start = -3)
pat_9_oment_met_pooled$transcript <- gsub('>', ':', pat_9_oment_met_pooled$transcript)
pat_9_oment_met_pooled$transcript <- paste0(pat_9_oment_met_pooled$chrom, ':', pat_9_oment_met_pooled$position, ':', pat_9_oment_met_pooled$transcript)
pat_9_oment_met_pooled <- pat_9_oment_met_pooled[, c('transcript', 'tumor_var_freq')]
colnames(pat_9_oment_met_pooled) <- c('Transcript', 'pat_9_oment_met')

pat_9_ovary_met_pooled <- pat_9_ovary_met_pooled[, c('position', 'chrom', 'transcript', 'tumor_var_freq')]
pat_9_ovary_met_pooled$transcript <- str_sub(pat_9_ovary_met_pooled$transcript, start = -3)
pat_9_ovary_met_pooled$transcript <- gsub('>', ':', pat_9_ovary_met_pooled$transcript)
pat_9_ovary_met_pooled$transcript <- paste0(pat_9_ovary_met_pooled$chrom, ':', pat_9_ovary_met_pooled$position, ':', pat_9_ovary_met_pooled$transcript)
pat_9_ovary_met_pooled <- pat_9_ovary_met_pooled[, c('transcript', 'tumor_var_freq')]
colnames(pat_9_ovary_met_pooled) <- c('Transcript', 'pat_9_ovary_met')

pat_9_plasma_pooled <- pat_9_plasma_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
pat_9_plasma_pooled$Transcript <- str_sub(pat_9_plasma_pooled$Transcript, start = -3)
pat_9_plasma_pooled$Transcript <- gsub('>', ':', pat_9_plasma_pooled$Transcript)
pat_9_plasma_pooled$Transcript <- paste0(pat_9_plasma_pooled$chrom, ':', pat_9_plasma_pooled$position, ':', pat_9_plasma_pooled$Transcript)
pat_9_plasma_pooled <- pat_9_plasma_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_9_plasma_pooled) <- c('Transcript', 'pat_9_plasma')

# put them together
pat_9_muts_pooled <- merge(pat_9_lymph_met_pooled, pat_9_oment_met_pooled, by = 'Transcript', all = TRUE)
pat_9_muts_pooled <- merge(pat_9_muts_pooled, pat_9_ovary_met_pooled, by = 'Transcript', all = TRUE)
pat_9_muts_pooled <- merge(pat_9_muts_pooled, pat_9_plasma_pooled, by = 'Transcript', all = TRUE)
pat_9_muts_pooled <- pat_9_muts_pooled[order(-pat_9_muts_pooled$pat_9_plasma), ]

rownames(pat_9_muts_pooled) <- pat_9_muts_pooled$Transcript
row_muts <- rownames(pat_9_muts_pooled)
pat_9_muts_pooled <- pat_9_muts_pooled[, -1]
#pat_9_muts_pooled <- sapply(pat_9_muts_pooled, function(x) ifelse (is.na(x), 0, x))
pat_9_muts_pooled <- as.data.frame(pat_9_muts_pooled)
rownames(pat_9_muts_pooled) <- row_muts
#plot figure to replace in manuscript


# custom layout for plots
layout_matrix <- matrix(c(2,2,3,3,3,
                          1,1,1,1,1), nrow = 2, ncol = 5, byrow = TRUE)
layout(mat = layout_matrix, 
       heights = c(2,2), 
       widths = c(1,1,1,1,1))
layout.show(3)


## CREATE MASTER COLOR SCALE and FIGURE??----
all_plasma <- c(pat_10_muts_pooled$pat_10_plasma, pat_2_muts_pooled$pat_2_plasma, pat_9_muts_pooled$pat_9_plasma)
plasma_cols <- rep('#000000', length(all_plasma))
plasma_cols <- color.scale(all_plasma, extremes = c('lightpink', 'red'), na.color = '#ffffff')

all_tumor <- c(pat_10_muts_pooled$pat_10_liv_met_1, pat_10_muts_pooled$pat_10_liv_met_2a, pat_10_muts_pooled$pat_10_liv_met_5, 
               pat_2_muts_pooled$pat_2_liv_met_1, pat_9_muts_pooled$pat_9_lymph_met, pat_9_muts_pooled$pat_9_oment_met, 
               pat_9_muts_pooled$pat_9_ovary_met)
tumor_cols <- rep('#000000', length(all_tumor))
tumor_cols <- color.scale(all_tumor, extremes = c('lightblue', 'blue'), na.color = '#ffffff')

master_cell_cols <- c(plasma_cols, tumor_cols)


# heatmap pat 10
pat_10_plasma_cols <- plasma_cols[1:27]
pat_10_tumor_cols <- tumor_cols[1:81]
pat_10_cell_cols <- c(pat_10_tumor_cols, pat_10_plasma_cols)
pat_10_cell_cols <- matrix(pat_10_cell_cols, nrow = 27, byrow = FALSE)
pat_10_cell_cols <- t(pat_10_cell_cols)
pat_10_cell_cols <- pat_10_cell_cols[c(4, 1:3), ]

pat_10_pooled_t <- data.frame(t(pat_10_muts_pooled))
pat_10_pooled_t <- pat_10_pooled_t[c(4, 1:3), ]


# plot it
# extra space
par(mar=c(6,15,6,2.1))

color2D.matplot(pat_10_pooled_t, cellcolors=pat_10_cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
legval<-seq(min(all_plasma, na.rm = TRUE),max(all_plasma, na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(3,-0.8,11,-0.4,round(c(min(all_plasma, na.rm = TRUE), max(all_plasma, na.rm = TRUE)),1),rect.col=legcol, font = 2)
mtext('Plasma', side=1, line=2.4, at=3.7, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(all_tumor, na.rm = TRUE),max(all_tumor, na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(13,-0.8,21,-0.4,round(c(min(all_tumor, na.rm = TRUE), max(all_tumor, na.rm = TRUE)),1),rect.col=legcol, font = 2)
mtext('Tumor', side=1, line=2.4, at=13.65, cex = 1.1, font = 2)
# add centered common text
mtext('Mutation Allele Frequency', side = 1, line = 4.1, at = 12, cex = 1.1, font = 2)

# add NA legend
color.legend(22.5, -0.8, 23, -0.4, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=3.1, at=23.9, cex = 1.1, font = 2)
legend(x=22.57,y=-.4,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
mut_col_end <- str_sub(mut_col_labels, -3)
mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 1, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver\nMet 1', 'Liver\nMet 2', 'Liver\nMet 5')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.5, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_1)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_1))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_2a)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_2a))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_10_muts_pooled$pat_10_liv_met_5)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_10_muts_pooled$pat_10_liv_met_5))), 
       pch = 16)

# add patient ID
mtext('Patient 10', side=1, line=1.7, at=-2.5, font = 2)


par(mar=c(5.1,4.1,4.1,2.1))


# heatmap pat 2 (only 1 mutation)
pat_2_plasma_cols <- plasma_cols[28]
pat_2_tumor_cols <- tumor_cols[82]
pat_2_cell_cols <- c(pat_2_tumor_cols, pat_2_plasma_cols)
pat_2_cell_cols <- matrix(pat_2_cell_cols, nrow = 1, byrow = FALSE)
pat_2_cell_cols <- t(pat_2_cell_cols)
pat_2_cell_cols <- pat_2_cell_cols[c(2,1), ]

pat_2_pooled_t <- data.frame(t(pat_2_muts_pooled))
pat_2_pooled_t <- pat_2_pooled_t[c(2,1), ]


# plot it
# extra space
par(mar=c(14.5,71.4,6,2.1)) # use these for single plot
par(mar=c(10.5,31.5,10.5,15.1)) # use these for laid out plots
color2D.matplot(data.frame(pat_2_pooled_t), cellcolors=pat_2_cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
# legval<-seq(min(pat_10_muts_pooled[, 4], na.rm = TRUE),max(pat_10_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
# legcol<-color.scale(legval, extremes = c('white', 'red'))
# color.legend(0,-0.9,10,-0.5,round(c(min(pat_10_muts_pooled[, 4], na.rm = TRUE), max(pat_10_muts_pooled[, 4], na.rm = TRUE)),1),rect.col=legcol)
# mtext('Plasma Mutant Allele Frequency', side=1, line=1.5, at=12.8)
# 
# # add tumor legend
# legval<-seq(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE),max(pat_10_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
# legcol<-color.scale(legval, extremes = c('white', 'blue'))
# color.legend(0,-1.7,10,-1.3,round(c(min(pat_10_muts_pooled[, 1:3], na.rm = TRUE), max(pat_10_muts_pooled[, 1:3], na.rm = TRUE)),1),rect.col=legcol)
# mtext('Tumor Mutant Allele Frequency', side=1, line=3.7, at=12.7)
# 
# # add NA legend
# color.legend(15.5, -1.7, 15.9, -1.3, legend = '', rect.col = '#00ff00')
# mtext('Mutation Not Present', side=1, line=3.7, at=17.7)
# legend(x=15.35,y=-1.15,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_2_muts_pooled)
mut_col_end <- str_sub(mut_col_labels, -3)
mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = 0.4, labels = mut_col_labels, tick = FALSE, cex.axis = 1, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Liver\nMet 1')
axis(2, at = c(0.6, 1.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.5, las = 1, font = 2)

# add patient ID
mtext('Patient 2', side=1, line=6.2, at=-3.0, font = 2)


par(mar=c(5.1,4.1,4.1,2.1))



# heatmap pat 9
pat_9_plasma_cols <- plasma_cols[29:41]
pat_9_tumor_cols <- tumor_cols[83:121]
pat_9_cell_cols <- c(pat_9_tumor_cols, pat_9_plasma_cols)
pat_9_cell_cols <- matrix(pat_9_cell_cols, nrow = 13, byrow = FALSE)
pat_9_cell_cols <- t(pat_9_cell_cols)
pat_9_cell_cols <- pat_9_cell_cols[c(4, 1:3), ]

pat_9_pooled_t <- data.frame(t(pat_9_muts_pooled))
pat_9_pooled_t <- pat_9_pooled_t[c(4, 1:3), ]


# plot it
# extra space
par(mar=c(6,45.5,6,2.1)) # use these for single plot
par(mar=c(6,21.5,6,2.1)) # use these for laid out plots

color2D.matplot(pat_9_pooled_t, cellcolors=pat_9_cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
# legval<-seq(min(pat_9_muts_pooled[, 4], na.rm = TRUE),max(pat_9_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
# legcol<-color.scale(legval, extremes = c('pink', 'red'))
# color.legend(0,-0.9,6,-0.5,round(c(min(pat_9_muts_pooled[, 4], na.rm = TRUE), max(pat_9_muts_pooled[, 4], na.rm = TRUE)),1),rect.col=legcol)
# mtext('Plasma Mutant Allele Frequency', side=1, line=1.5, at=7.4)
# 
# # add tumor legend
# legval<-seq(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE),max(pat_9_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
# legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
# color.legend(0,-1.7,6,-1.3,round(c(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE), max(pat_9_muts_pooled[, 1:3], na.rm = TRUE)),1),rect.col=legcol)
# mtext('Tumor Mutant Allele Frequency', side=1, line=3.7, at=7.35)
# 
# # add NA legend
# color.legend(8.7, -1.7, 8.9, -1.3, legend = '', rect.col = '#00ff00')
# mtext('Mutation Not Present', side=1, line=3.7, at=9.85)
# legend(x=8.65,y=-1.15,legend='',pch=20,bty="n",xpd = NA)

#plot labels THESE ARE NOT RIGHT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mut_col_labels <- rownames(pat_9_muts_pooled)
mut_col_end <- str_sub(mut_col_labels, -3)
mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_9_pooled_t)) - 0.6, labels = mut_col_end, tick = FALSE, cex.axis = 1, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Lymph\nMet', 'Omental\nMet', 'Ovary\nMet')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.5, las = 1, font = 2)

#add points for NA values
# lymph points
points(x = which(is.na(pat_9_muts_pooled$pat_9_lymph_met)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_9_muts_pooled$pat_9_lymph_met))), 
       pch = 16)
# omental points
points(x = which(is.na(pat_9_muts_pooled$pat_9_oment_met)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_9_muts_pooled$pat_9_oment_met))), 
       pch = 16)
# ovary points
points(x = which(is.na(pat_9_muts_pooled$pat_9_ovary_met)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_9_muts_pooled$pat_9_ovary_met))), 
       pch = 16)

# add patient ID
mtext('Patient 9', side=1, line=1.7, at=-3.0, font = 2)

par(mar=c(5.1,4.1,4.1,2.1))


pat_9_3 <- c(pat_9_muts_pooled$pat_9_lymph_met, pat_9_muts_pooled$pat_9_oment_met, pat_9_muts_pooled$pat_9_ovary_met)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_9_muts_pooled)[1] * dim(pat_9_muts_pooled)[2])
# plasma reds
cell_cols[241:320] <- color.scale(pat_9_muts_pooled[, 4], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:240] <- color.scale(pat_9_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 80, byrow = FALSE)
pat_9_pooled_t <- data.frame(t(pat_9_muts_pooled))
pat_9_pooled_t <- pat_9_pooled_t[c(4, 1:3), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(4, 1:3), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(6,15.5,6,12.1))
color2D.matplot(pat_9_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
legval<-seq(min(pat_9_muts_pooled[, 4], na.rm = TRUE),max(pat_9_muts_pooled[, 4], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(0,-0.9,30,-0.5,round(c(min(pat_9_muts_pooled[, 4], na.rm = TRUE), max(pat_9_muts_pooled[, 4], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=2.4, at=2.2, cex = 1.1, font = 2)

# add tumor legend
legval<-seq(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE),max(pat_9_muts_pooled[, 1:3], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(35,-0.9,65,-0.5,round(c(min(pat_9_muts_pooled[, 1:3], na.rm = TRUE), max(pat_9_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=2.4, at=37, cex = 1.1, font = 2)
mtext('Mutant Allele Frequency', side = 1, line = 4.3, at = 32.5, cex = 1.1, font = 2)

# add NA legend
color.legend(68, -0.9, 70, -0.5, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=3.0, at=72.9, cex = 1.1, font = 2)
legend(x=68.15,y=-0.47,legend='',pch=16,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_10_muts_pooled)
# mut_col_end <- str_sub(mut_col_labels, -3)
# mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
# mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_10_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.9, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Lymph\nMet', 'Omental\nMet', 'Ovary\nMet')
axis(2, at = c(0.6, 1.6, 2.6, 3.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)

#add points for NA values
# liver 1 points
points(x = which(is.na(pat_9_muts_pooled$pat_9_lymph_met)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_9_muts_pooled$pat_9_lymph_met))), 
       pch = 16)
# liver 2 points
points(x = which(is.na(pat_9_muts_pooled$pat_9_oment_met)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_9_muts_pooled$pat_9_oment_met))), 
       pch = 16)
# liver 5 points
points(x = which(is.na(pat_9_muts_pooled$pat_9_ovary_met)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_9_muts_pooled$pat_9_ovary_met))), 
       pch = 16)



par(mar=c(5.1,4.1,4.1,2.1))







## patient EMA----
# this patient has 7 mets and plasma
pat_ema_liv_met_2 <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP3_snv_annot.xlsx', 1)
# 361 left

pat_ema_oment_met_1 <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP4_snv_annot.xlsx', 1)
# 350 left

pat_ema_left_kidney_met <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP5_snv_annot.xlsx', 1)
# 306 left

pat_ema_right_kidney_met <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP6_snv_annot.xlsx', 1)
# 252 left

pat_ema_liv_met_1 <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP7_snv_annot.xlsx', 1)
# 364 left

pat_ema_heart_met <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP8_snv_annot.xlsx', 1)
# 194 left

pat_ema_oment_met_2 <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-CP11_snv_annot.xlsx', 1)
# 214 left

pat_ema_plasma <- mutation_process_2('N:/D13100 BGI analysis/Pt EMA/CP1-VS-UP2_snv_annot (low AF).xlsx', 1)
# 322 left

# what do they have in common?
all_pat_ema_genes <- Reduce(intersect, list(pat_ema_liv_met_2$Transcript, pat_ema_oment_met_1$Transcript, pat_ema_left_kidney_met$Transcript, 
                                            pat_ema_right_kidney_met$Transcript, pat_ema_liv_met_1$Transcript, pat_ema_heart_met$Transcript, 
                                            pat_ema_oment_met_2$Transcript)) #NONE in common across


# pooled vs plasma
pat_ema_pooled <- unique(c(pat_ema_liv_met_2$Transcript, pat_ema_oment_met_1$Transcript, pat_ema_left_kidney_met$Transcript, 
                           pat_ema_right_kidney_met$Transcript, pat_ema_liv_met_1$Transcript, pat_ema_heart_met$Transcript, 
                           pat_ema_oment_met_2$Transcript)) #1875

pat_ema_pooled_plasma <- intersect(pat_ema_pooled, pat_ema_plasma$Transcript) #38

# subset to pooled plasma and take a look
pat_ema_liv_met_2_pooled <- pat_ema_liv_met_2[pat_ema_liv_met_2$Transcript %in% pat_ema_pooled_plasma, ]
#9
pat_ema_oment_met_1_pooled <- pat_ema_oment_met_1[pat_ema_oment_met_1$Transcript %in% pat_ema_pooled_plasma, ]
#11
pat_ema_left_kidney_met_pooled <- pat_ema_left_kidney_met[pat_ema_left_kidney_met$Transcript %in% pat_ema_pooled_plasma, ]
#16
pat_ema_right_kidney_met_pooled <- pat_ema_right_kidney_met[pat_ema_right_kidney_met$Transcript %in% pat_ema_pooled_plasma, ]
#4
pat_ema_liv_met_1_pooled <- pat_ema_liv_met_1[pat_ema_liv_met_1$Transcript %in% pat_ema_pooled_plasma, ]
#14
pat_ema_heart_met_pooled <- pat_ema_heart_met[pat_ema_heart_met$Transcript %in% pat_ema_pooled_plasma, ]
#9
pat_ema_oment_met_2_pooled <- pat_ema_oment_met_2[pat_ema_oment_met_2$Transcript %in% pat_ema_pooled_plasma, ]
#7
pat_ema_plasma_pooled <- pat_ema_plasma[pat_ema_plasma$Transcript %in% pat_ema_pooled_plasma, ]
#38

# subset these to tumor AF for new figure
#pat_ema_liv_met_2_pooled <- pat_ema_liv_met_2_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
#pat_ema_liv_met_2_pooled$Transcript <- str_sub(pat_ema_liv_met_2_pooled$Transcript, start = -3)
#pat_ema_liv_met_2_pooled$Transcript <- gsub('>', ':', pat_ema_liv_met_2_pooled$Transcript)
#pat_ema_liv_met_2_pooled$Transcript <- paste0(pat_ema_liv_met_2_pooled$chrom, ':', pat_ema_liv_met_2_pooled$position, ':', pat_ema_liv_met_2_pooled$Transcript)
pat_ema_liv_met_2_pooled <- pat_ema_liv_met_2_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_liv_met_2_pooled) <- c('Transcript', 'pat_ema_liv_met_2')

#pat_ema_oment_met_1_pooled <- pat_ema_oment_met_1_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
#pat_ema_oment_met_1_pooled$Transcript <- str_sub(pat_ema_oment_met_1_pooled$Transcript, start = -3)
#pat_ema_oment_met_1_pooled$Transcript <- gsub('>', ':', pat_ema_oment_met_1_pooled$Transcript)
#pat_ema_oment_met_1_pooled$Transcript <- paste0(pat_ema_oment_met_1_pooled$chrom, ':', pat_ema_oment_met_1_pooled$position, ':', pat_ema_oment_met_1_pooled$Transcript)
pat_ema_oment_met_1_pooled <- pat_ema_oment_met_1_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_oment_met_1_pooled) <- c('Transcript', 'pat_ema_oment_met_1')

# pat_ema_left_kidney_met_pooled <- pat_ema_left_kidney_met_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_left_kidney_met_pooled$Transcript <- str_sub(pat_ema_left_kidney_met_pooled$Transcript, start = -3)
# pat_ema_left_kidney_met_pooled$Transcript <- gsub('>', ':', pat_ema_left_kidney_met_pooled$Transcript)
# pat_ema_left_kidney_met_pooled$Transcript <- paste0(pat_ema_left_kidney_met_pooled$chrom, ':', pat_ema_left_kidney_met_pooled$position, ':', pat_ema_left_kidney_met_pooled$Transcript)
pat_ema_left_kidney_met_pooled <- pat_ema_left_kidney_met_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_left_kidney_met_pooled) <- c('Transcript', 'pat_ema_left_kidney_met')

# pat_ema_right_kidney_met_pooled <- pat_ema_right_kidney_met_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_right_kidney_met_pooled$Transcript <- str_sub(pat_ema_right_kidney_met_pooled$Transcript, start = -3)
# pat_ema_right_kidney_met_pooled$Transcript <- gsub('>', ':', pat_ema_right_kidney_met_pooled$Transcript)
# pat_ema_right_kidney_met_pooled$Transcript <- paste0(pat_ema_right_kidney_met_pooled$chrom, ':', pat_ema_right_kidney_met_pooled$position, ':', pat_ema_right_kidney_met_pooled$Transcript)
pat_ema_right_kidney_met_pooled <- pat_ema_right_kidney_met_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_right_kidney_met_pooled) <- c('Transcript', 'pat_ema_right_kidney_met')

# pat_ema_liv_met_1_pooled <- pat_ema_liv_met_1_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_liv_met_1_pooled$Transcript <- str_sub(pat_ema_liv_met_1_pooled$Transcript, start = -3)
# pat_ema_liv_met_1_pooled$Transcript <- gsub('>', ':', pat_ema_liv_met_1_pooled$Transcript)
# pat_ema_liv_met_1_pooled$Transcript <- paste0(pat_ema_liv_met_1_pooled$chrom, ':', pat_ema_liv_met_1_pooled$position, ':', pat_ema_liv_met_1_pooled$Transcript)
pat_ema_liv_met_1_pooled <- pat_ema_liv_met_1_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_liv_met_1_pooled) <- c('Transcript', 'pat_ema_liv_met_1')

# pat_ema_heart_met_pooled <- pat_ema_heart_met_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_heart_met_pooled$Transcript <- str_sub(pat_ema_heart_met_pooled$Transcript, start = -3)
# pat_ema_heart_met_pooled$Transcript <- gsub('>', ':', pat_ema_heart_met_pooled$Transcript)
# pat_ema_heart_met_pooled$Transcript <- paste0(pat_ema_heart_met_pooled$chrom, ':', pat_ema_heart_met_pooled$position, ':', pat_ema_heart_met_pooled$Transcript)
pat_ema_heart_met_pooled <- pat_ema_heart_met_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_heart_met_pooled) <- c('Transcript', 'pat_ema_heart_met')

# pat_ema_oment_met_2_pooled <- pat_ema_oment_met_2_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_oment_met_2_pooled$Transcript <- str_sub(pat_ema_oment_met_2_pooled$Transcript, start = -3)
# pat_ema_oment_met_2_pooled$Transcript <- gsub('>', ':', pat_ema_oment_met_2_pooled$Transcript)
# pat_ema_oment_met_2_pooled$Transcript <- paste0(pat_ema_oment_met_2_pooled$chrom, ':', pat_ema_oment_met_2_pooled$position, ':', pat_ema_oment_met_2_pooled$Transcript)
pat_ema_oment_met_2_pooled <- pat_ema_oment_met_2_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_oment_met_2_pooled) <- c('Transcript', 'pat_ema_oment_met_2')

# pat_ema_plasma_pooled <- pat_ema_plasma_pooled[, c('position', 'chrom', 'Transcript', 'tumor_var_freq')]
# pat_ema_plasma_pooled$Transcript <- str_sub(pat_ema_plasma_pooled$Transcript, start = -3)
# pat_ema_plasma_pooled$Transcript <- gsub('>', ':', pat_ema_plasma_pooled$Transcript)
# pat_ema_plasma_pooled$Transcript <- paste0(pat_ema_plasma_pooled$chrom, ':', pat_ema_plasma_pooled$position, ':', pat_ema_plasma_pooled$Transcript)
pat_ema_plasma_pooled <- pat_ema_plasma_pooled[, c('Transcript', 'tumor_var_freq')]
colnames(pat_ema_plasma_pooled) <- c('Transcript', 'pat_ema_plasma')


# put them together
pat_ema_muts_pooled <- merge(pat_ema_heart_met_pooled, pat_ema_left_kidney_met_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_right_kidney_met_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_liv_met_1_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_liv_met_2_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_oment_met_1_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_oment_met_2_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- merge(pat_ema_muts_pooled, pat_ema_plasma_pooled, by = 'Transcript', all = TRUE)
pat_ema_muts_pooled <- pat_ema_muts_pooled[order(-pat_ema_muts_pooled$pat_ema_plasma), ]

rownames(pat_ema_muts_pooled) <- pat_ema_muts_pooled$Transcript
row_muts <- rownames(pat_ema_muts_pooled)
pat_ema_muts_pooled <- pat_ema_muts_pooled[, -1]
#pat_ema_muts_pooled <- sapply(pat_ema_muts_pooled, function(x) ifelse (is.na(x), 0, x))
pat_ema_muts_pooled <- as.data.frame(pat_ema_muts_pooled)
rownames(pat_ema_muts_pooled) <- row_muts
#plot figure to replace in manuscript

# 2d heatmap
# subset out plasma col
pat_ema_3 <- c(pat_ema_muts_pooled$pat_ema_heart_met, pat_ema_muts_pooled$pat_ema_left_kidney_met, pat_ema_muts_pooled$pat_ema_right_kidney_met, 
               pat_ema_muts_pooled$pat_ema_liv_met_1, pat_ema_muts_pooled$pat_ema_liv_met_2, pat_ema_muts_pooled$pat_ema_oment_met_1, 
               pat_ema_muts_pooled$pat_ema_oment_met_2)

#set colors, plasma has its own reds, pooled tumors blues
cell_cols<-rep("#000000",dim(pat_ema_muts_pooled)[1] * dim(pat_ema_muts_pooled)[2])
# plasma reds
cell_cols[1163:1328] <- color.scale(pat_ema_muts_pooled[, 8], extremes = c('lightpink', 'red'), na.color = '#ffffff')
# tumor blues
cell_cols[1:1162] <- color.scale(pat_ema_3, extremes = c('lightblue', 'blue'), na.color = '#ffffff')
cell_cols <- matrix(cell_cols, nrow = 166, byrow = FALSE)
pat_ema_pooled_t <- data.frame(t(pat_ema_muts_pooled))
pat_ema_pooled_t <- pat_ema_pooled_t[c(8, 1:7), ]

cell_cols <- t(cell_cols)
cell_cols <- cell_cols[c(8, 1:7), ]
# plot it
# extra space
par(mar=c(6,5.5,6,2.1))
#par(mar=c(5,10.5,9,2.1))

color2D.matplot(pat_ema_pooled_t, cellcolors=cell_cols, xlab = '', ylab = '', border='black', axes = FALSE)

# add plasma legend
legval<-seq(min(pat_ema_muts_pooled[, 8], na.rm = TRUE),max(pat_ema_muts_pooled[, 8], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightpink', 'red'))
color.legend(7,-1.4,67,-0.9,round(c(min(pat_ema_muts_pooled[, 8], na.rm = TRUE), max(pat_ema_muts_pooled[, 8], na.rm = TRUE)),2),rect.col=legcol)
mtext('Plasma', side=1, line=1.95, at=11.4, font = 2, cex = 1.1)

# add tumor legend
legval<-seq(min(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),max(pat_ema_muts_pooled[, 1:7], na.rm = TRUE),length.out = 100)
legcol<-color.scale(legval, extremes = c('lightblue', 'blue'))
color.legend(73,-1.4,133,-0.9,round(c(min(pat_ema_muts_pooled[, 1:3], na.rm = TRUE), max(pat_ema_muts_pooled[, 1:3], na.rm = TRUE)),2),rect.col=legcol)
mtext('Tumor', side=1, line=1.95, at=76.9, font = 2, cex = 1.1)
mtext('Mutant Allele Frequency', side = 1, line = 3.3, at = 68.3, cex = 1.1, font = 2)

# add NA legend
color.legend(136, -1.4, 139, -0.9, legend = '', rect.col = '#ffffff')
mtext('Mutation\n     Not Present', side=1, line=2.5, at=145.3, font = 2, cex = 1.1)
legend(x=135.6,y=-0.69,legend='',pch=20,bty="n",xpd = NA)

#plot labels
mut_col_labels <- rownames(pat_ema_muts_pooled)
mut_col_end <- str_sub(mut_col_labels, -3)
mut_col_labels <- gsub('.{3}$', '', mut_col_labels)
mut_col_labels <- paste0(mut_col_labels, '\n', mut_col_end)
axis(3, at = (1:ncol(pat_ema_pooled_t)) - 0.6, labels = mut_col_labels, tick = FALSE, cex.axis = 0.7, las = 2, font = 2)

mut_row_labels <- c('Plasma', 'Heart\nMet', 'Kidney\nMet L', 'Kidney\nMet R', 'Liver\nMet 1', 'Liver\nMet 2', 'Omental\nMet 1', 'Omental\nMet 2')
axis(2, at = c(0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6, 7.5), labels = rev(mut_row_labels), tick = FALSE, cex.axis = 1.1, las = 1, font = 2)


#add points for NA values
# heart points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_heart_met)) - 0.5, 
       y = rep(6.5, sum(is.na(pat_ema_muts_pooled$pat_ema_heart_met))), 
       pch = 20)
# left kidney points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_left_kidney_met)) - 0.5, 
       y = rep(5.5, sum(is.na(pat_ema_muts_pooled$pat_ema_left_kidney_met))), 
       pch = 20)
# right kidney points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_right_kidney_met)) - 0.5, 
       y = rep(4.5, sum(is.na(pat_ema_muts_pooled$pat_ema_right_kidney_met))), 
       pch = 20)
# liver met 1 points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_liv_met_1)) - 0.5, 
       y = rep(3.5, sum(is.na(pat_ema_muts_pooled$pat_ema_liv_met_1))), 
       pch = 20)
# liver met 2 points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_liv_met_2)) - 0.5, 
       y = rep(2.5, sum(is.na(pat_ema_muts_pooled$pat_ema_liv_met_2))), 
       pch = 20)
# oment met 1 points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_oment_met_1)) - 0.5, 
       y = rep(1.5, sum(is.na(pat_ema_muts_pooled$pat_ema_oment_met_1))), 
       pch = 20)
# oment met 2 points
points(x = which(is.na(pat_ema_muts_pooled$pat_ema_oment_met_2)) - 0.5, 
       y = rep(0.5, sum(is.na(pat_ema_muts_pooled$pat_ema_oment_met_2))), 
       pch = 20)







par(mar=c(5.1,4.1,4.1,2.1))









##create upset plots for sets before subsetting----
# patient 10
pat_10_liv_met_1_muts <- pat_10_liv_met_1$`Transcript change`
pat_10_liv_met_1_muts <- as.data.frame(pat_10_liv_met_1_muts)
pat_10_liv_met_1_muts$mutation <- rep(1, nrow(pat_10_liv_met_1_muts))
colnames(pat_10_liv_met_1_muts) <- c('mutation', 'Liver_Met_1')

pat_10_liv_met_2a_muts <- pat_10_liv_met_2a$`Transcript change`
pat_10_liv_met_2a_muts <- as.data.frame(pat_10_liv_met_2a_muts)
pat_10_liv_met_2a_muts$mutation <- rep(1, nrow(pat_10_liv_met_2a_muts))
colnames(pat_10_liv_met_2a_muts) <- c('mutation', 'Liver_Met_2a')

pat_10_liv_met_5_muts <- pat_10_liv_met_5$Transcript
pat_10_liv_met_5_muts <- as.data.frame(pat_10_liv_met_5_muts)
pat_10_liv_met_5_muts$mutation <- rep(1, nrow(pat_10_liv_met_5_muts))
colnames(pat_10_liv_met_5_muts) <- c('mutation', 'Liver_Met_5')

pat_10_plasma_muts <- pat_10_plasma$Transcript
pat_10_plasma_muts <- as.data.frame(pat_10_plasma_muts)
pat_10_plasma_muts$mutation <- rep(1, nrow(pat_10_plasma_muts))
colnames(pat_10_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_10_liv_met_1_muts, pat_10_liv_met_2a_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_10_liv_met_5_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_10_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 12, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(6)
purple_pal <- pal_material('purple', n= 8, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(2)
bar_colors <- c(orange_pal[1:6], purple_pal[1:2], 'darkgreen')


par(mar=c(5.1,4.1,4.1,2.1))

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Liver_Met_5 = 'white', Liver_Met_2a = 'white', 
                                                                 Liver_Met_1 = 'white')))), 
      intersections = list(list('Liver_Met_5', 'Plasma'), 
                           list('Liver_Met_2a', 'Plasma'), 
                           list('Liver_Met_1', 'Plasma'), 
                           list('Liver_Met_1', 'Liver_Met_2a'), 
                           list('Liver_Met_2a', 'Liver_Met_5'), 
                           list('Liver_Met_1', 'Liver_Met_5'), 
                           list('Liver_Met_1', 'Liver_Met_5', 'Plasma'), 
                           list('Liver_Met_1', 'Liver_Met_2a', 'Plasma'), 
                           list(colnames(upset_df))),
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray80', rep('firebrick4', 3)), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))


pat_10_liv_met_1_muts <- pat_10_liv_met_1$`Transcript change`
pat_10_liv_met_1_muts <- as.data.frame(pat_10_liv_met_1_muts)
pat_10_liv_met_1_muts$mutation <- rep(1, nrow(pat_10_liv_met_1_muts))
colnames(pat_10_liv_met_1_muts) <- c('mutation', 'Liver_Met_1')

pat_10_liv_met_2a_muts <- pat_10_liv_met_2a$`Transcript change`
pat_10_liv_met_2a_muts <- as.data.frame(pat_10_liv_met_2a_muts)
pat_10_liv_met_2a_muts$mutation <- rep(1, nrow(pat_10_liv_met_2a_muts))
colnames(pat_10_liv_met_2a_muts) <- c('mutation', 'Liver_Met_2a')

pat_10_liv_met_5_muts <- pat_10_liv_met_5$Transcript
pat_10_liv_met_5_muts <- as.data.frame(pat_10_liv_met_5_muts)
pat_10_liv_met_5_muts$mutation <- rep(1, nrow(pat_10_liv_met_5_muts))
colnames(pat_10_liv_met_5_muts) <- c('mutation', 'Liver_Met_5')

pat_10_plasma_muts <- pat_10_plasma$Transcript
pat_10_plasma_muts <- as.data.frame(pat_10_plasma_muts)
pat_10_plasma_muts$mutation <- rep(1, nrow(pat_10_plasma_muts))
colnames(pat_10_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_10_liv_met_1_muts, pat_10_liv_met_2a_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_10_liv_met_5_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_10_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
upset_df$mutation <- as.character(upset_df$mutation)

pat_10_pl_5 <- upset_df[(upset_df$Liver_Met_5 == 1 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 0 & upset_df$Liver_Met_2a == 0), ]
pat_10_pl_5 <- pat_10_pl_5$mutation
pat_10_liv_met_5_af <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_pl_5, ]
pat_10_liv_met_5_af <- pat_10_liv_met_5_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_pl_5, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_pl_5_af <- c(pat_10_liv_met_5_af, pat_10_plasma_af)
pat_10_pl_5_sd <- sd(pat_10_pl_5_af)
pat_10_pl_5_af <- median(pat_10_pl_5_af)

pat_10_pl_2a <- upset_df[(upset_df$Liver_Met_5 == 0 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 0 & upset_df$Liver_Met_2a == 1), ]
pat_10_pl_2a <- pat_10_pl_2a$mutation
pat_10_liv_met_2a_af <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_pl_2a, ]
pat_10_liv_met_2a_af <- pat_10_liv_met_2a_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_pl_2a, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_pl_2a_af <- c(pat_10_liv_met_2a_af, pat_10_plasma_af)
pat_10_pl_2a_sd <- sd(pat_10_pl_2a_af)
pat_10_pl_2a_af <- median(pat_10_pl_2a_af)

pat_10_pl_1 <- upset_df[(upset_df$Liver_Met_5 == 0 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 0), ]
pat_10_pl_1 <- pat_10_pl_1$mutation
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_pl_1, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_pl_1, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_pl_1_af <- c(pat_10_liv_met_1_af, pat_10_plasma_af)
pat_10_pl_1_sd <- sd(pat_10_pl_1_af)
pat_10_pl_1_af <- median(pat_10_pl_1_af)

pat_10_1_2a <- upset_df[(upset_df$Liver_Met_5 == 0 & upset_df$Plasma == 0 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 1), ]
pat_10_1_2a <- pat_10_1_2a$mutation
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_1_2a, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_liv_met_2a_af <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_1_2a, ]
pat_10_liv_met_2a_af <- pat_10_liv_met_2a_af$tumor_var_freq
pat_10_1_2a_af <- c(pat_10_liv_met_1_af, pat_10_liv_met_2a_af)
pat_10_1_2a_sd <- sd(pat_10_1_2a_af)
pat_10_1_2a_af <- median(pat_10_1_2a_af)

pat_10_2a_5 <- upset_df[(upset_df$Liver_Met_5 == 1 & upset_df$Plasma == 0 & upset_df$Liver_Met_1 == 0 & upset_df$Liver_Met_2a == 1), ]
pat_10_2a_5 <- pat_10_2a_5$mutation
pat_10_liv_met_5_af <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_2a_5, ]
pat_10_liv_met_5_af <- pat_10_liv_met_5_af$tumor_var_freq
pat_10_liv_met_2a_af <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_2a_5, ]
pat_10_liv_met_2a_af <- pat_10_liv_met_2a_af$tumor_var_freq
pat_10_2a_5_af <- c(pat_10_liv_met_5_af, pat_10_liv_met_2a_af)
pat_10_2a_5_sd <- sd(pat_10_2a_5_af)
pat_10_2a_5_af <- median(pat_10_2a_5_af)

pat_10_1_5 <- upset_df[(upset_df$Liver_Met_5 == 1 & upset_df$Plasma == 0 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 0), ]
pat_10_1_5 <- pat_10_1_5$mutation
pat_10_liv_met_5_af <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_1_5, ]
pat_10_liv_met_5_af <- pat_10_liv_met_5_af$tumor_var_freq
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_1_5, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_1_5_af <- c(pat_10_liv_met_5_af, pat_10_liv_met_1_af)
pat_10_1_5_sd <- sd(pat_10_1_5_af)
pat_10_1_5_af <- median(pat_10_1_5_af)

pat_10_1_5_pl <- upset_df[(upset_df$Liver_Met_5 == 1 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 0), ]
pat_10_1_5_pl <- pat_10_1_5_pl$mutation
pat_10_liv_met_5_af <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_1_5_pl, ]
pat_10_liv_met_5_af <- pat_10_liv_met_5_af$tumor_var_freq
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_1_5_pl, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_1_5_pl, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_1_5_pl_af <- c(pat_10_liv_met_5_af, pat_10_liv_met_1_af, pat_10_plasma_af)
pat_10_1_5_pl_sd <- sd(pat_10_1_5_pl_af)
pat_10_1_5_pl_af <- median(pat_10_1_5_pl_af)

pat_10_1_2a_pl <- upset_df[(upset_df$Liver_Met_5 == 0 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 1), ]
pat_10_1_2a_pl <- pat_10_1_2a_pl$mutation
pat_10_liv_met_2a_af <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_1_2a_pl, ]
pat_10_liv_met_2a_af <- pat_10_liv_met_2a_af$tumor_var_freq
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_1_2a_pl, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_1_2a_pl, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_1_2a_pl_af <- c(pat_10_liv_met_2a_af, pat_10_liv_met_1_af, pat_10_plasma_af)
pat_10_1_2a_pl_sd <- sd(pat_10_1_2a_pl_af)
pat_10_1_2a_pl_af <- median(pat_10_1_2a_pl_af)

pat_10_all <- upset_df[(upset_df$Liver_Met_5 == 1 & upset_df$Plasma == 1 & upset_df$Liver_Met_1 == 1 & upset_df$Liver_Met_2a == 1), ]
pat_10_all <- pat_10_all$mutation
pat_10_liv_met_2a_af <- pat_10_liv_met_2a[pat_10_liv_met_2a$`Transcript change` %in% pat_10_all, ]
pat_10_liv_met_2a_af <- pat_10_liv_met_2a_af$tumor_var_freq
pat_10_liv_met_1_af <- pat_10_liv_met_1[pat_10_liv_met_1$`Transcript change` %in% pat_10_all, ]
pat_10_liv_met_1_af <- pat_10_liv_met_1_af$tumor_var_freq
pat_10_plasma_af <- pat_10_plasma[pat_10_plasma$Transcript %in% pat_10_all, ]
pat_10_plasma_af <- pat_10_plasma_af$tumor_var_freq
pat_10_liv_met_5_af <- pat_10_liv_met_5[pat_10_liv_met_5$Transcript %in% pat_10_all, ]
pat_10_liv_met_5_af <- pat_10_liv_met_5_af$tumor_var_freq
pat_10_all_af <- c(pat_10_liv_met_2a_af, pat_10_liv_met_1_af, pat_10_plasma_af, pat_10_liv_met_5_af)
pat_10_all_sd <- sd(pat_10_all_af)
pat_10_all_af <- median(pat_10_all_af)


pat_10_afs <- c(pat_10_pl_5_af, pat_10_pl_2a_af, pat_10_pl_1_af, pat_10_1_2a_af, pat_10_2a_5_af, pat_10_1_5_af, 
                pat_10_1_5_pl_af, pat_10_1_2a_pl_af, pat_10_all_af)
pat_10_sds <- c(pat_10_pl_5_sd, pat_10_pl_2a_sd, pat_10_pl_1_sd, pat_10_1_2a_sd, pat_10_2a_5_sd, pat_10_1_5_sd, 
                pat_10_1_5_pl_sd, pat_10_1_2a_pl_sd, pat_10_all_sd)
pat_10_labels <- c('Liver Met 5\n+ Plasma', 'Liver Met 2a\n+ Plasma', 'Liver Met 1\n+ Plasma', 'Liver Met 1\n+ Liver Met 2a', 
                   'Liver Met 2a\n+ Liver Met 5', 'Liver Met 1\n+ Liver Met 5', 'Liver Met 1\n+ Liver Met 5\n+ Plasma', 
                   'Liver Met 1\n+ Liver Met 2a\n+ Plasma', 'All Samples')
pat_10_afs <- data.frame(pat_10_labels, pat_10_afs, pat_10_sds, bar_colors)
p<-ggplot(data=pat_10_afs, aes(x=pat_10_labels, y=pat_10_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_10_afs$pat_10_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors[c(7,8,1:4,5,6,9)]) +
  geom_errorbar(aes(ymin=pat_10_afs, ymax=pat_10_afs+pat_10_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=14, face = 'bold'),
                                                   axis.title=element_text(size=16,face="bold"))

p



# patient 2
pat_2_breast_2_muts <- pat_2_breast_2$`Transcript change`
pat_2_breast_2_muts <- as.data.frame(pat_2_breast_2_muts)
pat_2_breast_2_muts$mutation <- rep(1, nrow(pat_2_breast_2_muts))
colnames(pat_2_breast_2_muts) <- c('mutation', 'Breast_Primary')

pat_2_liv_met_1_muts <- pat_2_liv_met_1$`Transcript change`
pat_2_liv_met_1_muts <- as.data.frame(pat_2_liv_met_1_muts)
pat_2_liv_met_1_muts$mutation <- rep(1, nrow(pat_2_liv_met_1_muts))
colnames(pat_2_liv_met_1_muts) <- c('mutation', 'Liver_Met_1')

pat_2_liv_met_2_muts <- pat_2_liv_met_2$`Transcript change`
pat_2_liv_met_2_muts <- as.data.frame(pat_2_liv_met_2_muts)
pat_2_liv_met_2_muts$mutation <- rep(1, nrow(pat_2_liv_met_2_muts))
colnames(pat_2_liv_met_2_muts) <- c('mutation', 'Liver_Met_2')


pat_2_plasma_muts <- pat_2_plasma$`Transcript change`
pat_2_plasma_muts <- as.data.frame(pat_2_plasma_muts)
pat_2_plasma_muts$mutation <- rep(1, nrow(pat_2_plasma_muts))
colnames(pat_2_plasma_muts) <- c('mutation', 'Plasma')


#merge them all together
upset_df <- merge(pat_2_breast_2_muts, pat_2_liv_met_1_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_2_liv_met_2_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_2_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 6, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(4)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(1)
bar_colors <- c(blue_pal[1:4], orange_pal[1:4], purple_pal[1])


pat_2_upset <- upset(upset_df, set.metadata = list(data = metadata, 
                                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                                      colors = c(Plasma = 'gray50', Liver_Met_2 = 'white', Liver_Met_1 = 'white', 
                                                                                 Breast_Primary = 'white')))), 
                      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = TRUE, sets.x.label = 'Number of Mutations', sets.bar.color = c('gray60', rep('firebrick4', 2), 'deeppink'), 
                      matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
                      set_size.show = TRUE, set_size.numbers_size = 5, text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 1.3))


# patient 8
pat_8_axillary_met_muts <- pat_8_axillary$transcript
pat_8_axillary_met_muts <- as.data.frame(pat_8_axillary_met_muts)
pat_8_axillary_met_muts$mutation <- rep(1, nrow(pat_8_axillary_met_muts))
colnames(pat_8_axillary_met_muts) <- c('mutation', 'Axillary_Met')

pat_8_breast_1_muts <- pat_8_breast_1$transcript
pat_8_breast_1_muts <- as.data.frame(pat_8_breast_1_muts)
pat_8_breast_1_muts$mutation <- rep(1, nrow(pat_8_breast_1_muts))
colnames(pat_8_breast_1_muts) <- c('mutation', 'Breast_Primary_1')

pat_8_breast_2_met_muts <- pat_8_breast_2$transcript
pat_8_breast_2_met_muts <- as.data.frame(pat_8_breast_2_met_muts)
pat_8_breast_2_met_muts$mutation <- rep(1, nrow(pat_8_breast_2_met_muts))
colnames(pat_8_breast_2_met_muts) <- c('mutation', 'Breast_Primary_2')

pat_8_plasma_muts <- pat_8_plasma$transcript
pat_8_plasma_muts <- as.data.frame(pat_8_plasma_muts)
pat_8_plasma_muts$mutation <- rep(1, nrow(pat_8_plasma_muts))
colnames(pat_8_plasma_muts) <- c('mutation', 'Plasma')



#merge them all together
upset_df <- merge(pat_8_axillary_met_muts, pat_8_breast_1_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_8_breast_2_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_8_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 6, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(3)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(1)
bar_colors <- c(blue_pal[1:4], orange_pal[1:3], purple_pal[1])


pat_8_upset <- upset(upset_df, set.metadata = list(data = metadata, 
                                                   plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                                     colors = c(Plasma = 'gray50', Breast_Primary_2 = 'white', Breast_Primary_1 = 'white', 
                                                                                Axillary_Met = 'white')))), 
                     nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = TRUE, sets.x.label = 'Number of Mutations', sets.bar.color = c('gray60', rep('deeppink', 2), 'goldenrod4'), 
                     matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
                     set_size.show = TRUE, set_size.numbers_size = 5, text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 1.3))


# patient 9
pat_9_lymph_met_muts <- pat_9_lymph_met$transcript
pat_9_lymph_met_muts <- as.data.frame(pat_9_lymph_met_muts)
pat_9_lymph_met_muts$mutation <- rep(1, nrow(pat_9_lymph_met_muts))
colnames(pat_9_lymph_met_muts) <- c('mutation', 'Lymph_Met')

pat_9_oment_met_muts <- pat_9_oment_met$transcript
pat_9_oment_met_muts <- as.data.frame(pat_9_oment_met_muts)
pat_9_oment_met_muts$mutation <- rep(1, nrow(pat_9_oment_met_muts))
colnames(pat_9_oment_met_muts) <- c('mutation', 'Omental_Met')

pat_9_ovary_met_muts <- pat_9_ovary_met$transcript
pat_9_ovary_met_muts <- as.data.frame(pat_9_ovary_met_muts)
pat_9_ovary_met_muts$mutation <- rep(1, nrow(pat_9_ovary_met_muts))
colnames(pat_9_ovary_met_muts) <- c('mutation', 'Ovary_Met')

pat_9_plasma_muts <- pat_9_plasma$Transcript
pat_9_plasma_muts <- as.data.frame(pat_9_plasma_muts)
pat_9_plasma_muts$mutation <- rep(1, nrow(pat_9_plasma_muts))
colnames(pat_9_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_9_lymph_met_muts, pat_9_oment_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_ovary_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

rownames(upset_df) <- upset_df$mutation
row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:4])
randomnumber <- round(runif(4, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 6, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(4)
orange_pal <- pal_material('orange', n = 8, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(4)
purple_pal <- pal_material('purple', n= 4, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(2)
bar_colors <- c(orange_pal[1:4], purple_pal[1], 'darkgreen')

upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                                      colors = c(Plasma = 'gray60', Liver_Met_5 = 'white', Liver_Met_2a = 'white', 
                                                                 Liver_Met_1 = 'white')))), 
      intersections = list(list('Ovary_Met', 'Plasma'), 
                           list('Lymph_Met', 'Plasma'), 
                           list('Omental_Met', 'Plasma'), 
                           list('Omental_Met', 'Ovary_Met'), 
                           list('Omental_Met', 'Ovary_Met', 'Plasma'), 
                           list(colnames(upset_df))),
      nsets = 4, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', 
      sets.bar.color = c('gray60', 'goldenrod4', 'aquamarine3', 'chocolate3'), matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, 
      main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
      text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 2.0))

pat_9_oment_met_muts <- pat_9_oment_met$transcript
pat_9_oment_met_muts <- as.data.frame(pat_9_oment_met_muts)
pat_9_oment_met_muts$mutation <- rep(1, nrow(pat_9_oment_met_muts))
colnames(pat_9_oment_met_muts) <- c('mutation', 'Omental_Met')

pat_9_ovary_met_muts <- pat_9_ovary_met$transcript
pat_9_ovary_met_muts <- as.data.frame(pat_9_ovary_met_muts)
pat_9_ovary_met_muts$mutation <- rep(1, nrow(pat_9_ovary_met_muts))
colnames(pat_9_ovary_met_muts) <- c('mutation', 'Ovary_Met')

pat_9_lymph_met_muts <- pat_9_lymph_met$transcript
pat_9_lymph_met_muts <- as.data.frame(pat_9_lymph_met_muts)
pat_9_lymph_met_muts$mutation <- rep(1, nrow(pat_9_lymph_met_muts))
colnames(pat_9_lymph_met_muts) <- c('mutation', 'Lymph_Met')

pat_9_plasma_muts <- pat_9_plasma$Transcript
pat_9_plasma_muts <- as.data.frame(pat_9_plasma_muts)
pat_9_plasma_muts$mutation <- rep(1, nrow(pat_9_plasma_muts))
colnames(pat_9_plasma_muts) <- c('mutation', 'Plasma')

#merge them all together
upset_df <- merge(pat_9_oment_met_muts, pat_9_ovary_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_lymph_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_9_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
upset_df$mutation <- as.character(upset_df$mutation)

pat_9_pl_ov <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov <- pat_9_pl_ov$mutation
pat_9_ovary_met_af <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% pat_9_pl_ov, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$tumor_var_freq
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_pl_ov, ]
pat_9_plasma_af <- pat_9_plasma_af$tumor_var_freq
pat_9_pl_ov_af <- c(pat_9_ovary_met_af, pat_9_plasma_af)
pat_9_pl_ov_sd <- sd(pat_9_pl_ov_af)
pat_9_pl_ov_af <- median(pat_9_pl_ov_af)

pat_9_pl_ly <- upset_df[(upset_df$Omental_Met == 0 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 1), ]
pat_9_pl_ly <- pat_9_pl_ly$mutation
pat_9_lymph_met_af <- pat_9_lymph_met[pat_9_lymph_met$transcript %in% pat_9_pl_ly, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$tumor_var_freq
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_pl_ly, ]
pat_9_plasma_af <- pat_9_plasma_af$tumor_var_freq
pat_9_pl_ly_af <- c(pat_9_lymph_met_af, pat_9_plasma_af)
pat_9_pl_ly_sd <- sd(pat_9_pl_ly_af)
pat_9_pl_ly_af <- median(pat_9_pl_ly_af)


pat_9_pl_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 0 & upset_df$Lymph_Met == 0), ]
pat_9_pl_om <- pat_9_pl_om$mutation
pat_9_oment_met_af <- pat_9_oment_met[pat_9_oment_met$transcript %in% pat_9_pl_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$tumor_var_freq
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_pl_om, ]
pat_9_plasma_af <- pat_9_plasma_af$tumor_var_freq
pat_9_pl_om_af <- c(pat_9_oment_met_af, pat_9_plasma_af)
pat_9_pl_om_sd <- sd(pat_9_pl_om_af)
pat_9_pl_om_af <- median(pat_9_pl_om_af)

pat_9_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 0 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_ov_om <- pat_9_ov_om$mutation
pat_9_oment_met_af <- pat_9_oment_met[pat_9_oment_met$transcript %in% pat_9_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$tumor_var_freq
pat_9_ovary_met_af <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% pat_9_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$tumor_var_freq
pat_9_ov_om_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af)
pat_9_ov_om_sd <- sd(pat_9_ov_om_af)
pat_9_ov_om_af <- median(pat_9_ov_om_af)

pat_9_pl_ov_om <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 0), ]
pat_9_pl_ov_om <- pat_9_pl_ov_om$mutation
pat_9_ovary_met_af <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% pat_9_pl_ov_om, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$tumor_var_freq
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_pl_ov_om, ]
pat_9_plasma_af <- pat_9_plasma_af$tumor_var_freq
pat_9_oment_met_af <- pat_9_oment_met[pat_9_oment_met$transcript %in% pat_9_pl_ov_om, ]
pat_9_oment_met_af <- pat_9_oment_met_af$tumor_var_freq
pat_9_pl_ov_om_af <- c(pat_9_ovary_met_af, pat_9_plasma_af, pat_9_oment_met_af)
pat_9_pl_ov_om_sd <- sd(pat_9_pl_ov_om_af)
pat_9_pl_ov_om_af <- median(pat_9_pl_ov_om_af)

pat_9_all <- upset_df[(upset_df$Omental_Met == 1 & upset_df$Plasma == 1 & upset_df$Ovary_Met == 1 & upset_df$Lymph_Met == 1), ]
pat_9_all <- pat_9_all$mutation
pat_9_oment_met_af <- pat_9_oment_met[pat_9_oment_met$transcript %in% pat_9_all, ]
pat_9_oment_met_af <- pat_9_oment_met_af$tumor_var_freq
pat_9_ovary_met_af <- pat_9_ovary_met[pat_9_ovary_met$transcript %in% pat_9_all, ]
pat_9_ovary_met_af <- pat_9_ovary_met_af$tumor_var_freq
pat_9_plasma_af <- pat_9_plasma[pat_9_plasma$Transcript %in% pat_9_all, ]
pat_9_plasma_af <- pat_9_plasma_af$tumor_var_freq
pat_9_lymph_met_af <- pat_9_lymph_met[pat_9_lymph_met$transcript %in% pat_9_all, ]
pat_9_lymph_met_af <- pat_9_lymph_met_af$tumor_var_freq
pat_9_all_af <- c(pat_9_oment_met_af, pat_9_ovary_met_af, pat_9_plasma_af, pat_9_lymph_met_af)
pat_9_all_sd <- sd(pat_9_all_af)
pat_9_all_af <- median(pat_9_all_af)


pat_9_afs <- c(pat_9_pl_ov_af, pat_9_pl_ly_af, pat_9_pl_om_af, pat_9_ov_om_af, pat_9_pl_ov_om_af, pat_9_all_af)
pat_9_sds <- c(pat_9_pl_ov_sd, pat_9_pl_ly_sd, pat_9_pl_om_sd, pat_9_ov_om_sd, pat_9_pl_ov_om_sd, pat_9_all_sd)
pat_9_labels <- c('Ovary Met\n+ Plasma', 'Lymph Met\n+ Plasma', 'Omental Met\n+ Plasma', 'Omental Met\n+ Ovary Met', 
                  'Omental Met\n+ Ovary Met\n+ Lymph Met', 'All Samples')
pat_9_afs <- data.frame(pat_9_labels, pat_9_afs, pat_9_sds, bar_colors)
p<-ggplot(data=pat_9_afs, aes(x=pat_9_labels, y=pat_9_afs, fill = bar_colors)) +
  geom_bar(stat="identity", width = 0.5, color = 'black') + scale_x_discrete(limits=pat_9_afs$pat_9_labels) + ylab('Median Mutant Allele Frequency') +
  xlab('Intersections') + theme_bw() + 
  theme(panel.border = element_blank(), axis.line.x = element_line(), axis.line.y = element_line()) + 
  scale_fill_manual("legend", values = bar_colors[c(5,1:3,4,6)]) +
  geom_errorbar(aes(ymin=pat_9_afs, ymax=pat_9_afs+pat_9_sds), width=.2,
                position=position_dodge(.9)) + theme(axis.text=element_text(size=14, face = 'bold'),
                                                     axis.title=element_text(size=16,face="bold"))

p


# patient ema
pat_ema_heart_met_muts <- pat_ema_heart_met$Transcript
pat_ema_heart_met_muts <- as.data.frame(pat_ema_heart_met_muts)
pat_ema_heart_met_muts$mutation <- rep(1, nrow(pat_ema_heart_met_muts))
colnames(pat_ema_heart_met_muts) <- c('mutation', 'Heart_Met')

pat_ema_left_kidney_met_muts <- pat_ema_left_kidney_met$Transcript
pat_ema_left_kidney_met_muts <- as.data.frame(pat_ema_left_kidney_met_muts)
pat_ema_left_kidney_met_muts$mutation <- rep(1, nrow(pat_ema_left_kidney_met_muts))
colnames(pat_ema_left_kidney_met_muts) <- c('mutation', 'Left_Kidney_Met')

pat_ema_liv_met_1_muts <- pat_ema_liv_met_1$Transcript
pat_ema_liv_met_1_muts <- as.data.frame(pat_ema_liv_met_1_muts)
pat_ema_liv_met_1_muts$mutation <- rep(1, nrow(pat_ema_liv_met_1_muts))
colnames(pat_ema_liv_met_1_muts) <- c('mutation', 'Liver_Met_1')

pat_ema_liv_met_2_muts <- pat_ema_liv_met_2$Transcript
pat_ema_liv_met_2_muts <- as.data.frame(pat_ema_liv_met_2_muts)
pat_ema_liv_met_2_muts$mutation <- rep(1, nrow(pat_ema_liv_met_2_muts))
colnames(pat_ema_liv_met_2_muts) <- c('mutation', 'Liver_Met_2')

pat_ema_oment_met_1_muts <- pat_ema_oment_met_1$Transcript
pat_ema_oment_met_1_muts <- as.data.frame(pat_ema_oment_met_1_muts)
pat_ema_oment_met_1_muts$mutation <- rep(1, nrow(pat_ema_oment_met_1_muts))
colnames(pat_ema_oment_met_1_muts) <- c('mutation', 'Omental_Met_1')

pat_ema_oment_met_2_muts <- pat_ema_oment_met_2$Transcript
pat_ema_oment_met_2_muts <- as.data.frame(pat_ema_oment_met_2_muts)
pat_ema_oment_met_2_muts$mutation <- rep(1, nrow(pat_ema_oment_met_2_muts))
colnames(pat_ema_oment_met_2_muts) <- c('mutation', 'Omental_Met_2')

pat_ema_right_kidney_met_muts <- pat_ema_right_kidney_met$Transcript
pat_ema_right_kidney_met_muts <- as.data.frame(pat_ema_right_kidney_met_muts)
pat_ema_right_kidney_met_muts$mutation <- rep(1, nrow(pat_ema_right_kidney_met_muts))
colnames(pat_ema_right_kidney_met_muts) <- c('mutation', 'Right_Kidney_Met')

pat_ema_plasma_muts <- pat_ema_plasma$Transcript
pat_ema_plasma_muts <- as.data.frame(pat_ema_plasma_muts)
pat_ema_plasma_muts$mutation <- rep(1, nrow(pat_ema_plasma_muts))
colnames(pat_ema_plasma_muts) <- c('mutation', 'Plasma')


#merge them all together
upset_df <- merge(pat_ema_heart_met_muts, pat_ema_left_kidney_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_liv_met_1_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_liv_met_2_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_oment_met_1_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_oment_met_2_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_right_kidney_met_muts, by = 'mutation', all = TRUE)
upset_df <- merge(upset_df, pat_ema_plasma_muts, by = 'mutation', all = TRUE)
upset_df$mutation <- as.character(upset_df$mutation)

mut_idx <- is.na(upset_df$mutation)
table(mut_idx)
upset_df <- upset_df[!mut_idx, ]
rownames(upset_df) <- upset_df$mutation

row_muts <- rownames(upset_df)
upset_df <- upset_df[, -1]
upset_df <- sapply(upset_df, function(x) ifelse (is.na(x), 0, x))
upset_df <- as.data.frame(upset_df)
rownames(upset_df) <- row_muts


#set up dummy metadata for later features
sets_order <- colnames(upset_df[1:8])
randomnumber <- round(runif(8, min = 0, max = 90))
metadata <- as.data.frame(cbind(sets_order, randomnumber))
names(metadata) <- c("sets", "randomnumber")

blue_pal <- pal_material('blue', n = 10, alpha = 1, reverse = TRUE)
blue_pal <- blue_pal(8)
orange_pal <- pal_material('orange', n = 40, alpha = 1, reverse = TRUE)
orange_pal <- orange_pal(28)
purple_pal <- pal_material('purple', n= 30, alpha = 1, reverse = TRUE)
purple_pal <- purple_pal(18)
green_pal <- pal_material('green', n = 28, alpha = 1, reverse = TRUE)
green_pal <- green_pal(15)
bar_colors <- c(orange_pal[1:24], purple_pal[1:13])


upset(upset_df, set.metadata = list(data = metadata, 
                                    plots = list(list(type = 'matrix_rows', column = 'sets', 
                                    colors = c(Plasma = 'gray50', Heart_Met = 'white', Left_Kidney_Met = 'white', 
                                               Liver_Met_1 = 'white', Liver_Met_2 = 'white', 
                                               Omental_Met_1 = 'white', Omental_Met_2 = 'white', 
                                               Right_Kidney_Met = 'white')))), 
      intersections = list(list('Liver_Met_1', 'Plasma'), 
                           list('Heart_Met', 'Plasma'), 
                           list('Left_Kidney_Met', 'Plasma'), 
                           list('Liver_Met_2', 'Plasma'), 
                           list('Right_Kidney_Met', 'Plasma'), 
                           list('Omental_Met_2', 'Plasma'), 
                           list('Omental_Met_1', 'Plasma'), 
                           list('Omental_Met_1', 'Liver_Met_1'), 
                           list('Left_Kidney_Met', 'Omental_Met_1'), 
                           list('Liver_Met_2', 'Omental_Met_1'), 
                           list('Left_Kidney_Met', 'Liver_Met_1'), 
                           list('Left_Kidney_Met', 'Right_Kidney_Met'), 
                           list('Omental_Met_2', 'Omental_Met_1'), 
                           list('Right_Kidney_Met', 'Omental_Met_1'), 
                           list('Right_Kidney_Met', 'Liver_Met_1'), 
                           list('Liver_Met_2', 'Liver_Met_1'), 
                           list('Right_Kidney_Met', 'Liver_Met_2'), 
                           list('Left_Kidney_Met', 'Liver_Met_2'), 
                           list('Heart_Met', 'Liver_Met_2'), 
                           list('Heart_Met', 'Right_Kidney_Met'), 
                           list('Heart_Met', 'Left_Kidney_Met'), 
                           list('Omental_Met_2', 'Left_Kidney_Met'), 
                           list('Omental_Met_2', 'Liver_Met_2'), 
                           list('Omental_Met_2', 'Liver_Met_1'),
                           list('Liver_Met_2', 'Liver_Met_1', 'Plasma'), 
                           list('Left_Kidney_Met', 'Omental_Met_1', 'Plasma'), 
                           list('Left_Kidney_Met', 'Omental_Met_1', 'Liver_Met_2'), 
                           list('Heart_Met', 'Left_Kidney_Met', 'Plasma'), 
                           list('Left_Kidney_Met', 'Liver_Met_2', 'Plasma'), 
                           list('Omental_Met_2', 'Liver_Met_1', 'Plasma'), 
                           list('Right_Kidney_Met', 'Liver_Met_1', 'Plasma'), 
                           list('Heart_Met', 'Right_Kidney_Met', 'Left_Kidney_Met'), 
                           list('Heart_Met', 'Left_Kidney_Met', 'Omental_Met_1'), 
                           list('Omental_Met_2', 'Left_Kidney_Met', 'Omental_Met_2'), 
                           list('Right_Kidney_Met', 'Left_Kidney_Met', 'Liver_Met_1'), 
                           list('Omental_Met_2', 'Omental_Met_1', 'Liver_Met_1'), 
                           list('Left_Kidney_Met', 'Omental_Met_1', 'Liver_Met_1')),
      nsets = 8, nintersects = NA, sets = rev(sets_order), keep.order = FALSE, sets.x.label = 'Number of Mutations', sets.bar.color = c('gray60', rep('firebrick3', 2), 'chocolate3', rep('lawngreen', 2), 'chocolate3', 'red'), 
      matrix.color = 'midnightblue', matrix.dot.alpha = 0.8, main.bar.color = bar_colors, mainbar.y.label = 'Number of Mutations\nin Common', 
      set_size.show = TRUE, set_size.numbers_size = 5, text.scale = c(2.5, 1.5, 1.3, 1.3, 1.3, 1.3))



pdf(file = 'upsets.pdf', width = 11, height = 8.5, onefile = TRUE)
par(mfrow = c(5,1))
pat_10_upset
pat_2_upset
pat_8_upset
pat_9_upset
pat_ema_upset
dev.off()
par(mfrow = c(1,1))




abline(lm(pat_heart_plasma[, 2]~pat_heart_plasma[, 1]), col="red") # regression line (y~x) 
x <- rnorm(1000)
y <- rnorm(1000)
plot(x,y,xlim = c(0,1), ylim = c(0,1))
abline(lm(pat_ema_muts_pooled$pat_ema_plasma ~ pat_ema_muts_pooled$pat_ema_heart_met), col = 'red')
abline(lm(pat_ema_muts_pooled$pat_ema_plasma ~ pat_ema_muts_pooled$pat_ema_left_kidney_met), col = 'blue')




# red/black histograms
# patient 10
pat_10_liv_met_1_stats <- pat_10_liv_met_1[, c('Transcript change', 'tumor_var_freq')]
pat_10_liv_met_2a_stats <- pat_10_liv_met_2a[, c('Transcript change', 'tumor_var_freq')]
pat_10_liv_met_5_stats <- pat_10_liv_met_5[, c('Transcript', 'tumor_var_freq')]
colnames(pat_10_liv_met_5_stats)[1] <- 'Transcript change'
pat_10_all <- rbind(pat_10_liv_met_1_stats, pat_10_liv_met_2a_stats)
pat_10_all <- rbind(pat_10_all, pat_10_liv_met_5_stats)
pat_10_all <- pat_10_all[order(-pat_10_all$tumor_var_freq), ]
pat_10_all$color <- ifelse(pat_10_all$`Transcript change` %in% pat_10_plasma$Transcript, 'red', 'black')
barplot(pat_10_all$tumor_var_freq, col = pat_10_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 10\n(3Tumors)')

# patient 2
pat_2_breast_2_stats <- pat_2_breast_2[, c('Transcript change', 'tumor_var_freq')]
pat_2_liv_met_1_stats <- pat_2_liv_met_1[, c('Transcript change', 'tumor_var_freq')]
pat_2_liv_met_2_stats <- pat_2_liv_met_2[, c('Transcript change', 'tumor_var_freq')]
pat_2_all <- rbind(pat_2_breast_2_stats, pat_2_liv_met_1_stats)
pat_2_all <- rbind(pat_2_all, pat_2_liv_met_2_stats)
pat_2_all <- pat_2_all[order(-pat_2_all$tumor_var_freq), ]
pat_2_all$color <- ifelse(pat_2_all$`Transcript change` %in% pat_2_plasma$`Transcript change`, 'red', 'black')
barplot(pat_2_all$tumor_var_freq, col = pat_2_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 2\n(3 tumors)')

# patient 9
pat_9_lymph_met_stats <- pat_9_lymph_met[, c('transcript', 'tumor_var_freq')]
pat_9_oment_met_stats <- pat_9_oment_met[, c('transcript', 'tumor_var_freq')]
pat_9_ovary_met_stats <- pat_9_ovary_met[, c('transcript', 'tumor_var_freq')]
pat_9_all <- rbind(pat_9_lymph_met_stats, pat_9_oment_met_stats)
pat_9_all <- rbind(pat_9_all, pat_9_ovary_met_stats)
pat_9_all <- pat_9_all[order(-pat_9_all$tumor_var_freq), ]
pat_9_all$color <- ifelse(pat_9_all$transcript %in% pat_9_plasma$transcript, 'red', 'black')
barplot(pat_9_all$tumor_var_freq, col = pat_9_all$color, ylab = 'Mutant Allele Frequency in Tumor', 
        xlab = 'Variant', main = 'Patient 9\n(3 tumors)')

# patient ema
head(pat_ema_right_kidney_met)
pat_ema_heart_met_stats <- pat_ema_heart_met[, c('Transcript', 'tumor_var_freq')]
pat_ema_left_kidney_met_stats <- pat_ema_left_kidney_met[, c('Transcript', 'tumor_var_freq')]
pat_ema_liv_met_1_stats <- pat_ema_liv_met_1[, c('Transcript', 'tumor_var_freq')]
pat_ema_liv_met_2_stats <- pat_ema_liv_met_2[, c('Transcript', 'tumor_var_freq')]
pat_ema_oment_met_1_stats <- pat_ema_oment_met_1[, c('Transcript', 'tumor_var_freq')]
pat_ema_oment_met_2_stats <- pat_ema_oment_met_2[, c('Transcript', 'tumor_var_freq')]
pat_ema_right_kidney_met_stats <- pat_ema_right_kidney_met[, c('Transcript', 'tumor_var_freq')]
pat_ema_all <- rbind(pat_ema_heart_met_stats, pat_ema_left_kidney_met_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_liv_met_1_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_liv_met_2_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_oment_met_1_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_oment_met_2_stats)
pat_ema_all <- rbind(pat_ema_all, pat_ema_right_kidney_met_stats)
pat_ema_all <- pat_ema_all[order(-pat_ema_all$tumor_var_freq), ]
pat_ema_all$color <- ifelse(pat_ema_all$Transcript %in% pat_ema_plasma, 'red', 'black') # this doesn't work right now
barplot(pat_ema_all$tumor_var_freq, col = pat_ema_all$color)
table(pat_ema_all$color)

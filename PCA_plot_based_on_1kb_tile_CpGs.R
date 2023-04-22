# Drawing PCA plot for Nanopore methylation samples by tilling CpGs in 1kb length

# ------------------------------------------- Start -------------------------------------------

# Loading libraries

library(bsseq)
library(methylSig)
library(ggforce)
library(BiocParallel)
library(parallel)
library(data.table)

# --- Read modbam2bed outputs start ---

Chd1 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Chd1_CpGs.cpg.acc.bed")
colnames(Chd1) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                    "freq", "canon", "mod", "filt")

Ctr2 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ctr2_CpGs.cpg.acc.bed")
colnames(Ctr2) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                    "freq", "canon", "mod", "filt")

Ctr6_1 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ctr6_1_CpGs.cpg.acc.bed")
colnames(Ctr6_1) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Ctr6_2 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ctr6_2_CpGs.cpg.acc.bed")
colnames(Ctr6_2) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Ctr6_3 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ctr6_3_CpGs.cpg.acc.bed")
colnames(Ctr6_3) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Kmt2a_1 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Kmt2a_1_CpGs.cpg.acc.bed")
colnames(Kmt2a_1) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Kmt2a_2 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Kmt2a_2_CpGs.cpg.acc.bed")
colnames(Kmt2a_2) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Kmt2a_3 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Kmt2a_3_CpGs.cpg.acc.bed")
colnames(Kmt2a_3) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                   "freq", "canon", "mod", "filt")

Crebbp <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Crebbp_CpGs.cpg.acc.bed")
colnames(Crebbp) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                       "freq", "canon", "mod", "filt")

Ctr3 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ctr3_CpGs.cpg.acc.bed")
colnames(Ctr3) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                      "freq", "canon", "mod", "filt")

# Dnmt1 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Dnmt1_CpGs.cpg.acc.bed")
# colnames(Dnmt1) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                    "freq", "canon", "mod", "filt")

Ezh2 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Ezh2_CpGs.cpg.acc.bed")
colnames(Ezh2) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                     "freq", "canon", "mod", "filt")

Hdac6 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Hdac6_CpGs.cpg.acc.bed")
colnames(Hdac6) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                    "freq", "canon", "mod", "filt")

Kdm6a <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Kdm6a_CpGs.cpg.acc.bed")
colnames(Kdm6a) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                     "freq", "canon", "mod", "filt")

# --- Reading modbam2bed outputs end ---

# --- Generate BSseq objects start ---

# Why we assigned position argument to "start +1" ?, because this coordinate directly stands for actual position of cytosine in CpG. 
# Besides, when we do this, it makes difference in tile per number of CpGs when countOverlap() function 
# used for tiled and non-tiled BSseq objects.

Chd1$mod_and_canon <- Chd1$mod + Chd1$canon

Chd1_BSseq <- BSseq(chr = as.character(Chd1$chrom), pos = Chd1$start + 1, 
                    M = matrix(Chd1$mod), Cov = matrix(Chd1$mod_and_canon), sampleNames = c("Chd1"))

#

Ctr2$mod_and_canon <- Ctr2$mod + Ctr2$canon

Ctr2_BSseq <- BSseq(chr = as.character(Ctr2$chrom), pos = Ctr2$start + 1, 
                    M = matrix(Ctr2$mod), Cov = matrix(Ctr2$mod_and_canon), sampleNames = c("Ctr2"))

#

Ctr6_1$mod_and_canon <- Ctr6_1$mod + Ctr6_1$canon

Ctr6_1_BSseq <- BSseq(chr = as.character(Ctr6_1$chrom), pos = Ctr6_1$start +1, 
                   M = matrix(Ctr6_1$mod), Cov = matrix(Ctr6_1$mod_and_canon), sampleNames = c("Ctr6_1"))

#

Ctr6_2$mod_and_canon <- Ctr6_2$mod + Ctr6_2$canon

Ctr6_2_BSseq <- BSseq(chr = as.character(Ctr6_2$chrom), pos = Ctr6_2$start +1, 
                      M = matrix(Ctr6_2$mod), Cov = matrix(Ctr6_2$mod_and_canon), sampleNames = c("Ctr6_2"))

#

Ctr6_3$mod_and_canon <- Ctr6_3$mod + Ctr6_3$canon

Ctr6_3_BSseq <- BSseq(chr = as.character(Ctr6_3$chrom), pos = Ctr6_3$start +1, 
                      M = matrix(Ctr6_3$mod), Cov = matrix(Ctr6_3$mod_and_canon), sampleNames = c("Ctr6_3"))

#

Kmt2a_1$mod_and_canon <- Kmt2a_1$mod + Kmt2a_1$canon

Kmt2a_1_BSseq <- BSseq(chr = as.character(Kmt2a_1$chrom), pos = Kmt2a_1$start +1, 
                      M = matrix(Kmt2a_1$mod), Cov = matrix(Kmt2a_1$mod_and_canon), sampleNames = c("Kmt2a_1"))

#

Kmt2a_2$mod_and_canon <- Kmt2a_2$mod + Kmt2a_2$canon

Kmt2a_2_BSseq <- BSseq(chr = as.character(Kmt2a_2$chrom), pos = Kmt2a_2$start +1, 
                       M = matrix(Kmt2a_2$mod), Cov = matrix(Kmt2a_2$mod_and_canon), sampleNames = c("Kmt2a_2"))

#

Kmt2a_3$mod_and_canon <- Kmt2a_3$mod + Kmt2a_3$canon

Kmt2a_3_BSseq <- BSseq(chr = as.character(Kmt2a_3$chrom), pos = Kmt2a_3$start +1, 
                       M = matrix(Kmt2a_3$mod), Cov = matrix(Kmt2a_3$mod_and_canon), sampleNames = c("Kmt2a_3"))

#

Crebbp$mod_and_canon <- Crebbp$mod + Crebbp$canon

Crebbp_BSseq <- BSseq(chr = as.character(Crebbp$chrom), pos = Crebbp$start +1, 
                       M = matrix(Crebbp$mod), Cov = matrix(Crebbp$mod_and_canon), sampleNames = c("Crebbp"))

#

Ctr3$mod_and_canon <- Ctr3$mod + Ctr3$canon

Ctr3_BSseq <- BSseq(chr = as.character(Ctr3$chrom), pos = Ctr3$start +1, 
                      M = matrix(Ctr3$mod), Cov = matrix(Ctr3$mod_and_canon), sampleNames = c("Ctr3"))

#

# Dnmt1$mod_and_canon <- Dnmt1$mod + Dnmt1$canon

# Dnmt1_BSseq <- BSseq(chr = as.character(Dnmt1$chrom), pos = Dnmt1$start +1, 
                    # M = matrix(Dnmt1$mod), Cov = matrix(Dnmt1$mod_and_canon), sampleNames = c("Dnmt1"))

#

Ezh2$mod_and_canon <- Ezh2$mod + Ezh2$canon

Ezh2_BSseq <- BSseq(chr = as.character(Ezh2$chrom), pos = Ezh2$start +1, 
                     M = matrix(Ezh2$mod), Cov = matrix(Ezh2$mod_and_canon), sampleNames = c("Ezh2"))

#

Hdac6$mod_and_canon <- Hdac6$mod + Hdac6$canon

Hdac6_BSseq <- BSseq(chr = as.character(Hdac6$chrom), pos = Hdac6$start +1, 
                    M = matrix(Hdac6$mod), Cov = matrix(Hdac6$mod_and_canon), sampleNames = c("Hdac6"))

#

Kdm6a$mod_and_canon <- Kdm6a$mod + Kdm6a$canon

Kdm6a_BSseq <- BSseq(chr = as.character(Kdm6a$chrom), pos = Kdm6a$start +1, 
                     M = matrix(Kdm6a$mod), Cov = matrix(Kdm6a$mod_and_canon), sampleNames = c("Kdm6a"))

# --- Generate BSseq objects end ---

# --- Merging BSseq objects start ---

merged_BSseq <- combine(Chd1_BSseq, Ctr2_BSseq, Ctr6_1_BSseq,
                        Ctr6_2_BSseq, Ctr6_3_BSseq, Kmt2a_1_BSseq,
                        Kmt2a_2_BSseq, Kmt2a_3_BSseq, Crebbp_BSseq,
                        Ctr3_BSseq, Ezh2_BSseq,
                        Hdac6_BSseq, Kdm6a_BSseq)

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                                      "Kmt2a_3", "Crebbp", "Ctr3", "Ezh2", "Hdac6", "Kdm6a"))

nrow(merged_BSseq) # 2,314,224 shared 1kb tiles!

# --- Merging BSseq objects end ---

# --- Add sample name info to merged data start ---

metadata <- data.frame("Samples" = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                     "Kmt2a_3", "Crebbp", "Ctr3", "Ezh2", "Hdac6", "Kdm6a")) 
# sample description with metadata

pData(merged_BSseq)$metadata <- metadata
pData(merged_BSseq)

# --- Add sample name info to merged data end ---

# --- Filtering number of CpGs with 0 coverage in all samples start ---

merged_BSseq.cov <- getCoverage(merged_BSseq) # 21,860,443 CpGs in total
keep <- which(rowSums(merged_BSseq.cov) !=0) # 21,858,950 CpG has > 0 coverage in all samples
merged_BSseq_filtered <- merged_BSseq[keep,]
nrow(merged_BSseq_filtered) # 21,858,950 CpGs left!

# --- Filtering number of CpGs with 0 coverage in all samples end ---

# --- Tilling CpGs start ---

merged_BSseq_tiled <- tile_by_windows(bs = merged_BSseq_filtered, win_size = 1000) # 2,722,508 1kb length tiles!

# --- Tilling CpGs end ---

# --- Keeping tiles with at least 3 CpGs start ---

merged_BSseq_tiled_filtered <- merged_BSseq_tiled[rowSums(as.matrix(countOverlaps(merged_BSseq_tiled, merged_BSseq_filtered))) > 2,]

# 2,314,109 1kb lenght tiles left! 408,399 tiles are removed!

# --- Keeping tiles with at least 3 CpGs end ---

# --- Checking coverage of 1kb tile data start ---

## Number of tiles with 0 coverage in all samples

sum(rowSums(getCoverage(merged_BSseq_tiled_filtered)) == 0) # zero tile with zero coverage in all samples!
which(rowSums(getCoverage(merged_BSseq_tiled_filtered)) == 0)

# --- Checking coverage of 1kb tile data end ---

# --- Visualization of tiled CpG sites start ---

head(granges(merged_BSseq_tiled_filtered)) # see 1kb tiled coordinates
head(getCoverage(merged_BSseq_tiled_filtered)) # get coverage of data
round(colMeans(getCoverage(merged_BSseq_tiled_filtered)), 1) # average coverage of CpG tiles!

# --- Visualization of tiled CpG sites end ---

# --- get methylation matrix from filtered data with getMeth() function start ---

# Here raw data used instead of smoothed data

merged_BSseq_methylation <- getMeth(merged_BSseq_tiled_filtered, type = "raw")
nrow(merged_BSseq_methylation) # 2,314,109 tiles

# --- get methylation matrix from filtered data with getMeth() function end ---

# --- Filtering tiles with NA methylation values in at least one of the samples start ---

sum(rowSums(is.na(merged_BSseq_methylation)) > 0) # how many rows (tiles) have NA values
sum(apply(merged_BSseq_methylation, 1, anyNA)) # how many rows (tiles) have NA values
table(apply(merged_BSseq_methylation, 1, function(X) any(is.na(X)))) # how many rows (tiles) have NA values
which_NAs <- apply(merged_BSseq_methylation, 1, function(X) any(is.na(X))) # how many rows (tiles) have NA values
length(which(which_NAs)) # how many rows (tiles) have NA values
which(which_NAs) # Which rows (tiles) have NA values?
merged_BSseq_methylation[1238370,] # 1238370th tile has NA value
merged_BSseq_methylation[1238557,] # 1238557th tile has NA value

# There are 74,449 NA values, that is, there are 74,449 1kb tiles have NA methylation values in at least one of the samples

# Removing CpG tiles with NA values

merged_BSseq_methylation_filtered <- merged_BSseq_methylation[rowSums(is.na(merged_BSseq_methylation)) == 0,] # Removal NA values
sum(rowSums(is.na(merged_BSseq_methylation_filtered)) > 0)  # How many tiles have NA values?
nrow(merged_BSseq_methylation_filtered) # 2,239,660 tiles left across samples!

# --- Filtering tiles with NA methylation values in at least one of the samples end ---

# --- Sum of methylation values in each sample start ---

colSums(merged_BSseq_methylation_filtered)

# --- Sum of methylation values in each sample end ---

# --- PCA plot start ---

pca_data=prcomp(t(merged_BSseq_methylation_filtered))
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(merged_BSseq_methylation_filtered),
                         Samples = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                     "Kmt2a_3", "Crebbp", "Ctr3", "Ezh2", "Hdac6", "Kdm6a"))

ggplot(df_pca_data, aes(PC1,PC2, color = Samples,label=row.names(df_pca_data))) +
  geom_point(size=8)+ labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) +
  geom_text(nudge_x = 2.5,nudge_y = 2.5, size = 5) + theme_classic()

# --- PCA plot end ---
                   # ------------------------------------------- End -------------------------------------------

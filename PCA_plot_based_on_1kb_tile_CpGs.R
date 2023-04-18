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

Dnmt1 <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/Dnmt1_CpGs.cpg.acc.bed")
colnames(Dnmt1) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
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

Dnmt1$mod_and_canon <- Dnmt1$mod + Dnmt1$canon

Dnmt1_BSseq <- BSseq(chr = as.character(Dnmt1$chrom), pos = Dnmt1$start +1, 
                    M = matrix(Dnmt1$mod), Cov = matrix(Dnmt1$mod_and_canon), sampleNames = c("Dnmt1"))

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

# --- tilling each BSseq object to 1kb length size, so each tile has certain number of CpGs start ---

# tilling Chd1

Chd1_1kb_tiled <- tile_by_windows(bs = Chd1_BSseq, win_size = 1000)
nrow(Chd1_1kb_tiled) # 2,722,502 1kb-length tiles!

# Count overlaps between tiled genome and CpG level genome

table(countOverlaps(Chd1_1kb_tiled, Chd1_BSseq)) # hic CpG icermeyen 197,317 1kb tiled bölge var ama en fazla CpG iceren tiled ise
# 178 tane CpG iceriyor! bir tane CpG iceren ise 102,876 adet iken 3 tane CpG iceren tile sayısı ise 229,812' dir.
# Kasper' e not olsun diye start + 1 yapınca pozisyona gercekten de tile basına düsen CpG sayısı degisti; hic CpG icermeyen
# tile sayısı 197,360 oldu ama önceden 197,317 idi. En fazla CpG iceren tile ise 178 tane CpG iceriyordu simdi bu sayısı 177' ye
# düstü. 3 tane CpG iceren tile sayısı 229,812' den 229,911 tane oldu!

sum(table(countOverlaps(Chd1_1kb_tiled, Chd1_BSseq))) # this gives total number of 1kb tiles!

number_of_CpGs_per_tile <- as.data.frame(table(countOverlaps(Chd1_1kb_tiled, Chd1_BSseq)))
colnames(number_of_CpGs_per_tile) <- c("Number_of_CpGs", "Number_of_tiles_having_corresponding_number_of_CpGs")

ggplot(number_of_CpGs_per_tile, aes(x = Number_of_CpGs, y = Number_of_tiles_having_corresponding_number_of_CpGs)) +
  geom_bar(stat = "identity") +
  labs(x = "Number_of_CpGs", y = "Number_of_tiles_having_corresponding_number_of_CpGs", title = "1kb tile and CpG overlap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

number_of_CpGs_per_tile$number_of_total_CpGs <- (as.integer(number_of_CpGs_per_tile$Number_of_CpGs) - 1) * (as.integer(number_of_CpGs_per_tile$Number_of_tiles_having_corresponding_number_of_CpGs))
sum(number_of_CpGs_per_tile$number_of_total_CpGs) # 21,344,840 CpG var ama Chd1' de normalde 21,344,868 CpG var yani 28 CpG
# bir sekilde eksilmis oldu bunun sebebi bazen tile tam 1kb olmuyor 999bp de olabiliyor bu yüzden olabilir mi? Bazı CpG' leri
# bu yüzden kaybetmis olabiliriz! Burada bir sekilde tile sayısı ile tile' ların icerdigi CpG sayısını carparak toplam CpG sayısını
# yani original olarak BSseq olusturdugumda elde ettigim CpG sayısını elde etmeye calıstım fakat sanıyorum tile' ın mükemmel
# olarak 1kb yapılamamasından dolayı 28 tane CpG kaybettim!

# Keeping rows with at least 3 CpGs!, that is, keeping 1kb tiles with at least 3 CpGs

Chd1_1kb_tiled_filtered <- Chd1_1kb_tiled[rowSums(as.matrix(countOverlaps(Chd1_1kb_tiled, Chd1_BSseq))) > 2,] # 479,488 tane tile elenmesi
# lazım ve geriye 2,243,014 adet tile kalması gerekiyor..
nrow(Chd1_1kb_tiled_filtered) # 2,243,014 adet tile kaldı beklenildigi gibi..

table(countOverlaps(Chd1_1kb_tiled_filtered, Chd1_BSseq)) # Burada hic CpG iceremeyen tile olmaması gerekiyor!
# Gercekten de yok! Sanırım yaptıgım sey dogru!

# Visualization of number of CpGs and tiles after tile filtering having 2 CpGs and less!

number_of_CpGs_per_tile_2 <- as.data.frame(table(countOverlaps(Chd1_1kb_tiled_filtered, Chd1_BSseq)))
colnames(number_of_CpGs_per_tile_2) <- c("Number_of_CpGs", "Number_of_tiles_having_corresponding_number_of_CpGs")

ggplot(number_of_CpGs_per_tile_2, aes(x = Number_of_CpGs, y = Number_of_tiles_having_corresponding_number_of_CpGs)) +
  geom_bar(stat = "identity") +
  labs(x = "Number_of_CpGs", y = "Number_of_tiles_having_corresponding_number_of_CpGs", title = "1kb tile and CpG overlap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

head(granges(Chd1_1kb_tiled_filtered))

# tilling Ctr2 sample

Ctr2_1kb_tiled <- tile_by_windows(bs = Ctr2_BSseq, win_size = 1000)
nrow(Ctr2_1kb_tiled) # 2,722,508 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ctr2_1kb_tiled_filtered <- Ctr2_1kb_tiled[rowSums(as.matrix(countOverlaps(Ctr2_1kb_tiled, Ctr2_BSseq))) > 2,]
nrow(Ctr2_1kb_tiled_filtered) # 2,242,752 tiles left! 479,756 tiles removed!

# tilling Ctr6_1 sample

Ctr6_1_1kb_tiled <- tile_by_windows(bs = Ctr6_1_BSseq, win_size = 1000)
nrow(Ctr6_1_1kb_tiled) # 2,722,502 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ctr6_1_1kb_tiled_filtered <- Ctr6_1_1kb_tiled[rowSums(as.matrix(countOverlaps(Ctr6_1_1kb_tiled, Ctr6_1_BSseq))) > 2,]
nrow(Ctr6_1_1kb_tiled_filtered) # 2,242,969 tiles left!

# tilling Ctr6_2 sample

Ctr6_2_1kb_tiled <- tile_by_windows(bs = Ctr6_2_BSseq, win_size = 1000)
nrow(Ctr6_2_1kb_tiled) # 2,722,504 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ctr6_2_1kb_tiled_filtered <- Ctr6_2_1kb_tiled[rowSums(as.matrix(countOverlaps(Ctr6_2_1kb_tiled, Ctr6_2_BSseq))) > 2,]
nrow(Ctr6_2_1kb_tiled_filtered) # 2,243,425 tiles left!

# tilling Ctr6_3 sample

Ctr6_3_1kb_tiled <- tile_by_windows(bs = Ctr6_3_BSseq, win_size = 1000)
nrow(Ctr6_3_1kb_tiled) # 2,722,508 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ctr6_3_1kb_tiled_filtered <- Ctr6_3_1kb_tiled[rowSums(as.matrix(countOverlaps(Ctr6_3_1kb_tiled, Ctr6_3_BSseq))) > 2,]
nrow(Ctr6_3_1kb_tiled_filtered) # 2,311,551 tiles left!

# tilling Kmt2a_1 sample

Kmt2a_1_1kb_tiled <- tile_by_windows(bs = Kmt2a_1_BSseq, win_size = 1000)
nrow(Kmt2a_1_1kb_tiled) # 2,722,504 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Kmt2a_1_1kb_tiled_filtered <- Kmt2a_1_1kb_tiled[rowSums(as.matrix(countOverlaps(Kmt2a_1_1kb_tiled, Kmt2a_1_BSseq))) > 2,]
nrow(Kmt2a_1_1kb_tiled_filtered) # 2,243,282 tiles left!

# tilling Kmt2a_2 sample

Kmt2a_2_1kb_tiled <- tile_by_windows(bs = Kmt2a_2_BSseq, win_size = 1000)
nrow(Kmt2a_2_1kb_tiled) # 2,722,506 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Kmt2a_2_1kb_tiled_filtered <- Kmt2a_2_1kb_tiled[rowSums(as.matrix(countOverlaps(Kmt2a_2_1kb_tiled, Kmt2a_2_BSseq))) > 2,]
nrow(Kmt2a_2_1kb_tiled_filtered) # 2,722,506 tiles left!

# tilling Kmt2a_3 sample

Kmt2a_3_1kb_tiled <- tile_by_windows(bs = Kmt2a_3_BSseq, win_size = 1000)
nrow(Kmt2a_3_1kb_tiled) # 2,722,508 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Kmt2a_3_1kb_tiled_filtered <- Kmt2a_3_1kb_tiled[rowSums(as.matrix(countOverlaps(Kmt2a_3_1kb_tiled, Kmt2a_3_BSseq))) > 2,]
nrow(Kmt2a_3_1kb_tiled_filtered) # 2,312,208 tiles left!

# tilling Crebbp sample

Crebbp_1kb_tiled <- tile_by_windows(bs = Crebbp_BSseq, win_size = 1000)
nrow(Crebbp_1kb_tiled) # 2,722,502 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Crebbp_1kb_tiled_filtered <- Crebbp_1kb_tiled[rowSums(as.matrix(countOverlaps(Crebbp_1kb_tiled, Crebbp_BSseq))) > 2,]
nrow(Crebbp_1kb_tiled_filtered) # 2,243,613 tiles left!

# tilling Ctr3 sample

Ctr3_1kb_tiled <- tile_by_windows(bs = Ctr3_BSseq, win_size = 1000)
nrow(Ctr3_1kb_tiled) # 2,722,502 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ctr3_1kb_tiled_filtered <- Ctr3_1kb_tiled[rowSums(as.matrix(countOverlaps(Ctr3_1kb_tiled, Ctr3_BSseq))) > 2,]
nrow(Ctr3_1kb_tiled_filtered) # 2,243,458 tiles left!

# tilling Dnmt1 sample

Dnmt1_1kb_tiled <- tile_by_windows(bs = Dnmt1_BSseq, win_size = 1000)
nrow(Dnmt1_1kb_tiled) # 2,722,496 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Dnmt1_1kb_tiled_filtered <- Dnmt1_1kb_tiled[rowSums(as.matrix(countOverlaps(Dnmt1_1kb_tiled, Dnmt1_BSseq))) > 2,]
nrow(Dnmt1_1kb_tiled_filtered) # 2,243,022 tiles left!

# tilling Ezh2 sample

Ezh2_1kb_tiled <- tile_by_windows(bs = Ezh2_BSseq, win_size = 1000)
nrow(Ezh2_1kb_tiled) # 2,722,495 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Ezh2_1kb_tiled_filtered <- Ezh2_1kb_tiled[rowSums(as.matrix(countOverlaps(Ezh2_1kb_tiled, Ezh2_BSseq))) > 2,]
nrow(Ezh2_1kb_tiled_filtered) # 2,243,195 tiles left!

# tilling Hdac6 sample

Hdac6_1kb_tiled <- tile_by_windows(bs = Hdac6_BSseq, win_size = 1000)
nrow(Hdac6_1kb_tiled) # 2,722,502 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Hdac6_1kb_tiled_filtered <- Hdac6_1kb_tiled[rowSums(as.matrix(countOverlaps(Hdac6_1kb_tiled, Hdac6_BSseq))) > 2,]
nrow(Hdac6_1kb_tiled_filtered) # 2,243,302 tiles left!

# tilling Kdm6a sample

Kdm6a_1kb_tiled <- tile_by_windows(bs = Kdm6a_BSseq, win_size = 1000)
nrow(Kdm6a_1kb_tiled) # 2,722,493 tiles!

# Count overlaps between tiled genome and CpG level genome and keep tiles with at least 3 CpGs

Kdm6a_1kb_tiled_filtered <- Kdm6a_1kb_tiled[rowSums(as.matrix(countOverlaps(Kdm6a_1kb_tiled, Kdm6a_BSseq))) > 2,]
nrow(Kdm6a_1kb_tiled_filtered) # 2,243,316 tiles left!

# --- tile each BSseq object to 1kb length size, so each tile has certain number of CpGs end ---

# --- Merging BSseq object start ---

merged_BSseq <- combine(Chd1_1kb_tiled_filtered, Ctr2_1kb_tiled_filtered, Ctr6_1_1kb_tiled_filtered,
                        Ctr6_2_1kb_tiled_filtered, Ctr6_3_1kb_tiled_filtered, Kmt2a_1_1kb_tiled_filtered,
                        Kmt2a_2_1kb_tiled_filtered, Kmt2a_3_1kb_tiled_filtered, Crebbp_1kb_tiled_filtered,
                        Ctr3_1kb_tiled_filtered, Dnmt1_1kb_tiled_filtered, Ezh2_1kb_tiled_filtered,
                        Hdac6_1kb_tiled_filtered, Kdm6a_1kb_tiled_filtered)

# Collapsedata

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                                      "Kmt2a_3", "Crebbp", "Ctr3", "Dnmt1", "Ezh2", "Hdac6", "Kdm6a"))

nrow(merged_BSseq) # 2,314,224 shared 1kb tiles!

# Add sample name info to merged data

metadata <- data.frame("Samples" = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                     "Kmt2a_3", "Crebbp", "Ctr3", "Dnmt1", "Ezh2", "Hdac6", "Kdm6a")) 
# sample description with metadata

pData(merged_BSseq)$metadata <- metadata
pData(merged_BSseq)

# --- Merging BSseq object end ---

# --- Coverage filtering for 1kb tile data start ---

## Number of tiles with 0 coverage in all samples

sum(rowSums(getCoverage(merged_BSseq)) == 0) # 50 tane tile' ın coverage' ı sıfır! geriye 2,314,174 tile kalması lazım!

# Filtering tiles with 0 coverage in all samples

merged_BSseq.cov <- getCoverage(merged_BSseq)
keep <- which(rowSums(merged_BSseq.cov) !=0)
merged_BSseq_2 <- merged_BSseq[keep,]
nrow(merged_BSseq_2) # Beklenildigi gibi 2,314,174 adet tile kaldı geriye ve bunların hepsi bir sekilde coverage' a sahip'ler!

# --- Coverage filtering for 1kb tile data end ---

# --- Visualization of tiled CpG sites start ---

head(granges(merged_BSseq_2)) # see 1kb tiled coordinates
head(getCoverage(merged_BSseq_2)) # get coverage of data
round(colMeans(getCoverage(merged_BSseq_2)), 1) # average coverage of CpG tiles!

# --- Visualization of tiled CpG sites end ---

# --- get methylation matrix from filtered data with getMeth() function start ---

# Here raw data used instead of smoothed data

merged_BSseq_methylation <- getMeth(merged_BSseq_2, type = "raw")
nrow(merged_BSseq_methylation) # 2,314,174 tiles

# --- get methylation matrix from filtered data with getMeth() function end ---

# --- Filtering non-overlapping tiles between samples which returns NA values in rows start ---

sum(rowSums(is.na(merged_BSseq_methylation)) > 0) # how many rows (tiles) have NA values
sum(apply(merged_BSseq_methylation, 1, anyNA)) # how many rows (tiles) have NA values
table(apply(merged_BSseq_methylation, 1, function(X) any(is.na(X)))) # how many rows (tiles) have NA values
which_NAs <- apply(merged_BSseq_methylation, 1, function(X) any(is.na(X))) # how many rows (tiles) have NA values
length(which(which_NAs)) # how many rows (tiles) have NA values
which(which_NAs) # Which rows (tiles) have NA values?
merged_BSseq_methylation[2116940,] # 2116940th tile has NA value
merged_BSseq_methylation[2139948,] # 2139948th tile has NA value

# There are 74,667 NA values, that is, there are 74,667 non-overlapping CpG tiles

# Removing CpG tiles with NA values

merged_BSseq_methylation_filtered <- merged_BSseq_methylation[rowSums(is.na(merged_BSseq_methylation)) == 0,] # Removal NA values
sum(rowSums(is.na(merged_BSseq_methylation_filtered)) > 0)  # How many tiles have NA values?
nrow(merged_BSseq_methylation_filtered) # 2,239,507 overlapping tiles left across samples!

# --- Filtering non-overlapping tiles between samples which returns NA values in rows end ---

# --- PCA plot start ---

pca_data=prcomp(t(merged_BSseq_methylation_filtered))
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(merged_BSseq_methylation_filtered),
                         Samples = c("Chd1", "Ctr2", "Ctr6_1", "Ctr6_2", "Ctr6_3", "Kmt2a_1", "Kmt2a_2",
                                       "Kmt2a_3", "Crebbp", "Ctr3", "Dnmt1", "Ezh2", "Hdac6", "Kdm6a"))

ggplot(df_pca_data, aes(PC1,PC2, color = Samples,label=row.names(df_pca_data))) +
  geom_point(size=8)+ labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) +
  geom_text(nudge_x = 2.5,nudge_y = 2.5, size = 5) + theme_classic()

# --- PCA plot end ---

# ------------------------------------------- End -------------------------------------------
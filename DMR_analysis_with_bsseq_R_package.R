# --------------------- Kmt2a start ---------------------
library(rtracklayer)
library(bsseq)
library(BiocParallel)
library(parallel)
library(tidyr)
library(tidyselect)
library(broom)
library(bsseqData)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(scales)
library(colorspace)
library(dichromat)
library(BRGenomics)

# --- Reading modbam2bed files  start ---

Kmt2a_1_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Kmt2a_1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Kmt2a_1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")

nrow(Kmt2a_1_bed_file) # 21,267,038 CpG loci var.

Kmt2a_2_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Kmt2a_2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Kmt2a_2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")

nrow(Kmt2a_2_bed_file) # 21,263,598 CpG loci var

# Control samples

Ctr6_1_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Ctr6_1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Ctr6_1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                               "freq", "canon", "mod", "filt")

Ctr6_2_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Ctr6_2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Ctr6_2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                               "freq", "canon", "mod", "filt")

Ctr2_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Ctr2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Ctr2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")

Ctr3_bed_file <- read.table("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Ctr3_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
colnames(Ctr3_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")

# --- Reading modbam2bed files  end ---

# --- Coverage calculation for each CpG by sum canonical and modified CpGs and then generation of BSseq objects for each sample start ---

Kmt2a_1_bed_file$mod_and_canon <- Kmt2a_1_bed_file$mod + Kmt2a_1_bed_file$canon

# BSseq object for Kmt2a_1

Kmt2a_1_BSseq <- BSseq(chr = as.character(Kmt2a_1_bed_file$chrom), pos = Kmt2a_1_bed_file$start +1, 
                       M = matrix(Kmt2a_1_bed_file$mod), Cov = matrix(Kmt2a_1_bed_file$mod_and_canon), sampleNames = c("Kmt2a_1"))

sampleNames(Kmt2a_1_BSseq) <- "Kmt2a_1"

# Coverage calculation by counting canonical and modified CpGs

Kmt2a_2_bed_file$mod_and_canon <- Kmt2a_2_bed_file$mod + Kmt2a_2_bed_file$canon

# BSseq object for Kmt2a_2

Kmt2a_2_BSseq <- BSseq(chr = as.character(Kmt2a_2_bed_file$chrom), pos = Kmt2a_2_bed_file$start +1, 
                       M = matrix(Kmt2a_2_bed_file$mod), Cov = matrix(Kmt2a_2_bed_file$mod_and_canon), sampleNames = c("Kmt2a_2"))

sampleNames(Kmt2a_2_BSseq) <- "Kmt2a_2"

# BSseq object for Ctr6_1

Ctr6_1_bed_file$mod_and_canon <- Ctr6_1_bed_file$mod + Ctr6_1_bed_file$canon

Ctr6_1_BSseq <- BSseq(chr = as.character(Ctr6_1_bed_file$chrom), pos = Ctr6_1_bed_file$start+1, 
                      M = matrix(Ctr6_1_bed_file$mod), Cov = matrix(Ctr6_1_bed_file$mod_and_canon), sampleNames = c("Ctr6_1"))

sampleNames(Ctr6_1_BSseq) <- "Ctr6_1"

# BSseq object for Ctr6_2

Ctr6_2_bed_file$mod_and_canon <- Ctr6_2_bed_file$mod + Ctr6_2_bed_file$canon

# Ctr6_2 icin BSseq obj yaratılması

Ctr6_2_BSseq <- BSseq(chr = as.character(Ctr6_2_bed_file$chrom), pos = Ctr6_2_bed_file$start+1, 
                      M = matrix(Ctr6_2_bed_file$mod), Cov = matrix(Ctr6_2_bed_file$mod_and_canon), sampleNames = c("Ctr6_2"))

sampleNames(Ctr6_2_BSseq) <- "Ctr6_2"

# BSseq object for Ctr2

Ctr2_bed_file$mod_and_canon <- Ctr2_bed_file$mod + Ctr2_bed_file$canon

Ctr2_BSseq <- BSseq(chr = as.character(Ctr2_bed_file$chrom), pos = Ctr2_bed_file$start +1, 
                    M = matrix(Ctr2_bed_file$mod), Cov = matrix(Ctr2_bed_file$mod_and_canon), sampleNames = c("Ctr2"))

sampleNames(Ctr2_BSseq) <- "Ctr2"

# BSseq object for Ctr3

Ctr3_bed_file$mod_and_canon <- Ctr3_bed_file$mod + Ctr3_bed_file$canon

Ctr3_BSseq <- BSseq(chr = as.character(Ctr3_bed_file$chrom), pos = Ctr3_bed_file$start +1, 
                    M = matrix(Ctr3_bed_file$mod), Cov = matrix(Ctr3_bed_file$mod_and_canon), sampleNames = c("Ctr3"))

sampleNames(Ctr3_BSseq) <- "Ctr3"

# --- Coverage calculation for each CpG by sum canonical and modified CpGs and then generation of BSseq objects for each sample start ---

# --- Combining BSseq objects start ---

merged_BSseq <- combine(Ctr6_1_BSseq, Ctr6_2_BSseq, Ctr2_BSseq, Ctr3_BSseq, Kmt2a_1_BSseq, Kmt2a_2_BSseq)

# --- Combining BSseq objects end ---

# --- Collapsedata BSseq objects start ---

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("Control1", "Control2", "Control3", "Control4", "Kmt2a_1", "Kmt2a_2"))

# --- Collapsedata BSseq objects end ---

# --- Smooting data start ---

start.time <- Sys.time()

merged_BSseq_smooted <- BSmooth(
  BSseq = merged_BSseq, 
  BPPARAM = MulticoreParam(workers = 8), 
  verbose = TRUE)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# --- Smooting data end ---

# --- Some statistics on data start ---

head(getCoverage(merged_BSseq_smooted, type = "Cov"), n = 6) # Get coverage for samples

length(merged_BSseq_smooted) # How many CpGs in all samples

round(colMeans(getCoverage(merged_BSseq_smooted)), 1) # Average coverage in all samples

sum(rowSums(getCoverage(merged_BSseq_smooted)) == 0) # Number of CpGs with 0 coverage in all samples

# --- Some statistics on data end ---

# --- Computing t-statistics start ---

# Creating metadata for all samples

metadata_2 <- data.frame("Samples" = c("Control", "Control", "Control","Control",
                                       "Kmt2a", "Kmt2a"))

pData(merged_BSseq_smooted)$metadata_2 <- metadata_2
pData(merged_BSseq_smooted)

# remove CpGs with little or no coverage
# we will only keep CpGs where at least 1 Kmt2a samples and at least 2 Control samples have at least 10x in coverage

Kmt2a.cov <- getCoverage(merged_BSseq_smooted)
keepLoci.ex <- which(rowSums(Kmt2a.cov[, merged_BSseq_smooted$metadata_2 == "Kmt2a"] >= 10) >= 1 &
                       rowSums(Kmt2a.cov[, merged_BSseq_smooted$metadata_2 == "Control"] >= 10) >= 2)

length(keepLoci.ex) 

merged_BSseq_smooted_filtered <- merged_BSseq_smooted[keepLoci.ex,]

# Computing t-statistics

merged_BSseq_smooted_filtered.tstat <- BSmooth.tstat(merged_BSseq_smooted_filtered, 
                                                     group1 = c("Kmt2a_1", "Kmt2a_2"),
                                                     group2 = c("Ctr6_1", "Ctr6_2", "Ctr2", "Ctr3"),
                                                     estimate.var = "group2",
                                                     local.correct = TRUE,
                                                     verbose = TRUE)

# Kasper recommended that "k" argument in BSmooth.tstat() function should be set 21. This argument affects variance of DMRs across samples!.
# k = 21 is standard deviation smooting. When you changed this cutoff to k=21, variance in DMRs across samples should be changed! Check it
# in DMR plots for example MORC4 gene!

plot(merged_BSseq_smooted_filtered.tstat)

# --- Computing t-statistics end ---

# Once t-statistics have been computed, we can compute differentially methylated regions (DMRs) by thresholding the t-statistics.
# Cutoff is chosen by looking at the quantiles of the t-statistics (for the entire genome).

# --- Finding DMRs start ---

dmrs0 <- dmrFinder(merged_BSseq_smooted_filtered.tstat, qcutoff = c(0.025, 0.975)) # cutoff for quantiles of the t-statistics.
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

# Kasper recommended that qcutoff should be removed because when we set this cutoff we say that always we want to have five percent of
# genome should have DMRs, equally, that equal number of hypermethylation and hypomethylation! So this give us always DMRs. Because 0.025
# is left tail of DMRs histogram but 0.975 is right tail of histogram of DMRs. So this always guarantee us to have DMRs even if they are
# not true positives! So, we should not use this cutoff, instead we should use cutoff = c(-4.6, 4.6) in dmrFinder() function to find DMRs.
# This absolute cutoff guarantee us imbalance number between hypermethylated and hypomethylated DMRs! So we'll have different number of
# hypermethylation and hypomethylation.

# Here, we filter out DMRs that do not have at least 3 CpGs in them
# and at least a mean difference (across the DMR) in methylation between control and Kmt2a of at least 0.1.

# --- Finding DMRs end ---

nrow(dmrs) # 1075 DMRs!
table(dmrs$direction) # 597 hypermethylation and 478 hypomethylation

write.table(dmrs, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/Kmt2a_DMRs.txt", quote = F, row.names = F, sep =  "\t")

# --- Plotting by assigning some specific colors to each group start ---

pData <- pData(merged_BSseq_smooted_filtered)
pData$col <- c(rep("#7570b3", each=4), rep("#d95f02", each=2)) # Purple is control, Orange is Kmt2a
pData(merged_BSseq_smooted_filtered) <- pData
pData(merged_BSseq_smooted_filtered)

# --- Plotting by assigning some specific colors to each group end ---

# --- Plotting only one DMR start ---

plotRegion(merged_BSseq_smooted_filtered, dmrs[3,], extend = 5000, addRegions = dmrs) # Plotting DMR which is in 3rd row.

# --- Plotting only one DMR end ---

# --- Plotting all DMRs start ---

pdf(file = "Kmt2a_DMRs.pdf", width = 10, height = 5)
plotManyRegions(merged_BSseq_smooted_filtered, dmrs[1:1075,], extend = 5000, 
                addRegions = dmrs) # plotting 1075 DMRs
dev.off()

# --- Plotting all DMRs end ---

# --- Heatmap for DMRs start ---

# Get coordinates of all DMRs as GRanges object

DMR_regions_2 <- GRanges(seqnames = c(as.character(dmrs$chr)), 
                         ranges = IRanges(start = dmrs$start, end = dmrs$end))

# get smoothed methylation values of all DMRs

Smoothed_Methylation_level_of_DMR_regions_2 <- getMeth(merged_BSseq_smooted_filtered, DMR_regions_2, type = "smooth", what = "perRegion")
head(Smoothed_Methylation_level_of_DMR_regions_2)
range(Smoothed_Methylation_level_of_DMR_regions_2) # Lowest methylation value is 0.05602372, biggest one is = 0.97391033
write.csv(Smoothed_Methylation_level_of_DMR_regions_2, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/Kmt2a_Smooted_methylation_values_of_DMRs.csv", quote = F, row.names = F)

# Draw heatmap of all DMRs based on smoothed methylation values

pheatmap::pheatmap(Smoothed_Methylation_level_of_DMR_regions_2, color=colorRampPalette(c("navy", "gray90", "red"))(50), cellwidth = 30, cellheight = 0.45)

# --- Heatmap for DMRs end ---

# Finding EM genes with DMRs

EM_genes <- fread("/media/ko/New Volume/Downloads/The Epigenetic Machinery.csv")
Kmt2a_DMRs_and_nearest_gene_overlap <- read.delim("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/nearest_geneids.txt", header = F)
table(intersect(EM_genes$`Gene Name`, toupper(Kmt2a_DMRs_and_nearest_gene_overlap$V2)))

# "CHD8", "HDGF", "HR", "KDM5C", "MORC4", "PHF8", "PRDM5", "PWWP2A", "SCML2", "SND1", "TDRD3",
# "TET2", and "TRIM66" EM genes overlap DMRs or nearest DMRs!

# Visualization of particular DMR region

# --- Showing EM genes' DMRs with heatmap start ---

# 1. Subset DMRs of EM genes according to their start positions

subset_DMRs <- subset(dmrs, dmrs$start %in% c("52243603", "87908595", "70569594", "152233983", "139870471", "151499108",
                                              "65933001", "43683117", "161161019", "161162080", "161162598",
                                              "28927798", "87569312", "87912448", "133545392", "109494582"))

# 2. Generation of GRanges coordinates for DMRs

GRanges <- GRanges(seqnames = c(as.character(subset_DMRs$chr)),
                   ranges = IRanges(start = subset_DMRs$start, end = subset_DMRs$end))

# 3. Obtaining smoothed methylation values of DMRs

methylation_values <- getMeth(merged_BSseq_smooted_filtered, GRanges, type = "smooth", what = "perRegion")

# 4. Manually adding DMR names to methylation matrix

write.csv(methylation_values, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/smoothed_methylation_values_of_DMRs.csv", row.names = F, quote = F)

# 4.1. get methylation value of one DMR region as GRanges object

TET2_DMR_region_1 <- subset(dmrs, dmrs$start == "133545392")

TET2_DMR_region_3.1 <- GRanges(seqnames = c(as.character(TET2_DMR_region_1$chr)),
                           ranges = IRanges(start = TET2_DMR_region_1$start, end = TET2_DMR_region_1$end))

# 4.2. Get methylation value of this particular DMR and find methylation values in previous metyhlation matrix (methylation_values) to label DMRs!

getMeth(merged_BSseq_smooted_filtered, TET2_DMR_region_3.1, type = "smooth", what = "perRegion")

# 5. Uploading modified metylation matrix

modified_methylation_matrix <- read.csv("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/smoothed_methylation_values_of_DMRs.csv", row.names = 1)

# 5. Drawing heatmap based on methylation values

pheatmap::pheatmap(modified_methylation_matrix, color=colorRampPalette(c("navy", "gray90", "red"))(50), cellwidth = 50, cellheight = 30)

# --- Showing EM genes' DMRs with heatmap end ---

# --- plotting gene/promoter/enhancer tracks with plotRegion() function in R package bsseq start ---

# 1. Convert bed file containing tracks of enhancer, promoter etc. to a GRanges object

# MORC4 DMR coordinates: chrX:139870471-139871803

MORC4_tracks <- "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/EM_genes_tracks_bed_files/MORC4_tracks.bed"
MORC4_tracks <- import(MORC4_tracks)

# 2. Subsetting tracks

CpG_island_88_ <- subset(MORC4_tracks, MORC4_tracks$name == "CpG:_88")
EM10E0930002 <- subset(MORC4_tracks, MORC4_tracks$name == "EM10E0930002")
EM10E0930003 <- subset(MORC4_tracks, MORC4_tracks$name == "EM10E0930003")
EM10E0930004 <- subset(MORC4_tracks, MORC4_tracks$name == "EM10E0930004")
EM10E0930005 <- subset(MORC4_tracks, MORC4_tracks$name == "EM10E0930005")
OREG1827292 <- subset(MORC4_tracks, MORC4_tracks$name == "OREG1827292")
OREG0959289 <- subset(MORC4_tracks, MORC4_tracks$name == "OREG0959289")
OREG0592150 <- subset(MORC4_tracks, MORC4_tracks$name == "OREG0592150")
CACCCn <- subset(MORC4_tracks, MORC4_tracks$name == "(CACCC)n")
intron2 <- subset(MORC4_tracks, MORC4_tracks$name == "intron2")
intron1 <- subset(MORC4_tracks, MORC4_tracks$name == "intron1")
exon2 <- subset(MORC4_tracks, MORC4_tracks$name == "exon2")
exon1 <- subset(MORC4_tracks, MORC4_tracks$name == "exon1")

# 3. Merging all tracks into one GRanges object

MORC4_tracks_merged <- GRangesList(CpG_island_88 = CpG_island_88_, Enhancer = EM10E0930002, Enhancer = EM10E0930003,
                                   Enhancer = EM10E0930004, Promoter = EM10E0930005, Mtf2_binding_site = OREG1827292,
                                   Myod1_binding_site = OREG0959289, E2F3_binding_site = OREG0592150, CACCCn_repeat = CACCCn,
                                   Intron2 = intron2, Intron1 = intron1, Exon2 = exon2, Exon1 = exon1)
# 4. Get DMRs of Morc4 gene

MORC4_DMR_region <- subset(dmrs, dmrs$start == "139870471")

# 5. Visualization of Morc4 DMR region with tracks

plotRegion(BSseq = merged_BSseq_smooted_filtered, region = MORC4_DMR_region, extend = 5000, addRegions = MORC4_DMR_region, main = "MORC4", annoTrack = MORC4_tracks_merged, cex.anno = 1, lwd = 2.5, mainWithWidth = F, regionCol = "#1b9e77")

# in plotRegion() function, region argument is interval of genomic coordinates but addRegions argument is used for highlighting of DMRs on plot.
# annoTrack argument is used for adding promoter, enhancer, transcription factor binding sites as tracks on plot.

# When you have consequentially distributed DMR regions, that is, they are next to next each other in the same gene
# how to display two different DMRs of the same gene in the same plot?
# For instance, SCML2 gene has three different DMRs which are relatively close each other but bsseq package outputs each of those DMRs separately.
# When you need to display for all three DMRs in the sample plot, you can execute following command below:

# SCML2 DMR coordinates 1: chrX:161161019-161161113
# SCML2 DMR coordinates 2: chrX:161162080-161162192
# SCML2 DMR coordinates 3: chrX:161162598-161163226

# Get tracks for SCML2 gene

SCML2_tracks <- "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/EM_genes_tracks_bed_files/SCML2_tracks.bed"
SCML2_tracks <- import(SCML2_tracks)

# Subsetting each tracks to one

AT_rich <- subset(SCML2_tracks, SCML2_tracks$name == "AT_rich")
T_rich <- subset(SCML2_tracks, SCML2_tracks$name == "T-rich")
EM10E0931270 <- subset(SCML2_tracks, SCML2_tracks$name == "EM10E0931270")
EM10E0931271 <- subset(SCML2_tracks, SCML2_tracks$name == "EM10E0931271")
EM10E0931272 <- subset(SCML2_tracks, SCML2_tracks$name == "EM10E0931272")

# Generation of GRanges object that has all track info

SCML2_merged <- GRangesList(AT_rich_repeat = AT_rich, T_rich_repeat = T_rich,
                            Enhancer = EM10E0931270, Enhancer = EM10E0931271, Enhancer = EM10E0931272)

# Subsetting DMRs and generation of GRanges object of all three DMRs

SCML2_region_1 <- subset(dmrs, dmrs$start == "161161019")
SCML2_region_2 <- subset(dmrs, dmrs$start == "161162080")
SCML2_region_3 <- subset(dmrs, dmrs$start == "161162598")

SCML2_region_merged <- GRangesList(SCML2_region_1, SCML2_region_2, SCML2_region_3)
# SCML2_region_merged object has three different DMRs, separately.

# for region argument in plotRegion() function, generation of genomic interval for DMR of interest!

SCML2_region_merged_2 <- GRanges(seqnames = "chrX", ranges = IRanges(start = 161156019, end = 161168226))
# SCML2_region_merged_2 object contains interval of three different DMRs of SCML2

plotRegion(BSseq = merged_BSseq_smooted_filtered, region = SCML2_region_merged_2, addRegions = SCML2_region_merged, main = "SCML2", annoTrack = SCML2_merged, cex.anno = 1, lwd = 2.5, mainWithWidth = F, regionCol = "#1b9e77")

# Now, three different DMRs in the same gene are displayed in the same plot instead of displaying in the different plots.

# --------------------- Kmt2a end ---------------------

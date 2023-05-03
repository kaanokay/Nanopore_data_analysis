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

# Coverage calculation for each CpG by sum canonical and modified CpGs

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

# Combining BSseq objects

merged_BSseq <- combine(Ctr6_1_BSseq, Ctr6_2_BSseq, Ctr2_BSseq, Ctr3_BSseq, Kmt2a_1_BSseq, Kmt2a_2_BSseq)

# Collapsedata

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("Control1", "Control2", "Control3", "Control4", "Kmt2a_1", "Kmt2a_2"))

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

head(getCoverage(merged_BSseq_smooted, type = "Cov"), n = 6) # Get coverage for samples

length(merged_BSseq_smooted) # How many CpGs in all samples

round(colMeans(getCoverage(merged_BSseq_smooted)), 1) # Average coverage in all samples

sum(rowSums(getCoverage(merged_BSseq_smooted)) == 0) # ## Number of CpGs with 0 coverage in all samples

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

merged_BSseq_smooted_filtered.tstat <- BSmooth.tstat(merged_BSseq_smooted_filtered, 
                                                     group1 = c("Kmt2a_1", "Kmt2a_2"),
                                                     group2 = c("Ctr6_1", "Ctr6_2", "Ctr2", "Ctr3"),
                                                     estimate.var = "group2",
                                                     local.correct = TRUE,
                                                     verbose = TRUE)

plot(merged_BSseq_smooted_filtered.tstat)

# --- Computing t-statistics end ---

# Once t-statistics have been computed, we can compute differentially methylated regions (DMRs) by thresholding the t-statistics.
# Here we use a cutoff of 4.6, which was chosen by looking at the quantiles of the t-statistics (for the entire genome).

# Finding DMRs

dmrs0 <- dmrFinder(merged_BSseq_smooted_filtered.tstat, qcutoff = c(0.025, 0.975)) # cutoff for quantiles of the t-statistics.
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

# Here, we filter out DMRs that do not have at least 3 CpGs in them and
# at least a mean difference (across the DMR) in methylation between control and Kmt2a of at least 0.1.

nrow(dmrs) # 1075 DMRs!
head(dmrs, n = 3)
table(dmrs$direction) # 1075 DMRs! 597 hypermethylation + 478 hypomethylation

write.table(dmrs, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/Kmt2a_DMRs.txt", quote = F, row.names = F, sep =  "\t")

# Plotting by assigning some specific colors to each group

pData <- pData(merged_BSseq_smooted_filtered)
pData$col <- c(rep("#7570b3", each=4), rep("#d95f02", each=2)) # Purple is control, Orange is Kmt2a
pData(merged_BSseq_smooted_filtered) <- pData
pData(merged_BSseq_smooted_filtered)

# Plotting only one DMR start

plotRegion(merged_BSseq_smooted_filtered, dmrs[3,], extend = 5000, addRegions = dmrs) # Plotting DMR which is in 3rd row.

# Plotting only one DMR end

# Plotting all DMRs start

pdf(file = "Kmt2a_DMRs.pdf", width = 10, height = 5)
plotManyRegions(merged_BSseq_smooted_filtered, dmrs[1:1075,], extend = 5000, 
                addRegions = dmrs)
dev.off()

# Plotting all DMRs end

# Heatmap for Differentially methylated CpGs start

# Get coordinates of all DMRs as GRanges object

DMR_regions_2 <- GRanges(seqnames = c(as.character(dmrs$chr)), 
                         ranges = IRanges(start = dmrs$start, end = dmrs$end))

# get smoothed methylation values of DMRs start

Smoothed_Methylation_level_of_DMR_regions_2 <- getMeth(merged_BSseq_smooted_filtered, DMR_regions_2, type = "smooth", what = "perRegion")
head(Smoothed_Methylation_level_of_DMR_regions_2)
range(Smoothed_Methylation_level_of_DMR_regions_2) # Lowest methylation value is 0.05602372, biggest one is = 0.97391033
write.csv(Smoothed_Methylation_level_of_DMR_regions_2, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/Kmt2a_Smooted_methylation_values_of_DMRs.csv", quote = F, row.names = F)

# get smoothed methylation values of DMRs end

# Draw heatmap of all DMRs based on smoothed methylation values

pheatmap::pheatmap(Smoothed_Methylation_level_of_DMR_regions_2, color=colorRampPalette(c("navy", "gray90", "red"))(50), cellwidth = 30, cellheight = 0.45)

# --- Heatmap for Differentially methylated CpGs end ---

# Finding EM genes with DMRs

EM_genes <- fread("/media/ko/New Volume/Downloads/The Epigenetic Machinery.csv")
Kmt2a_DMRs_and_nearest_gene_overlap <- read.delim("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/nearest_geneids.txt", header = F)
table(intersect(EM_genes$`Gene Name`, toupper(Kmt2a_DMRs_and_nearest_gene_overlap$V2)))

# CHD8, HDGF, HR, KDM5C, MORC4, PHF8, PRDM5, PWWP2A, SCML2, SND1, TDRD3, TET2, TRIM66 EM genlerinde DMR var!

# Visualization of particular DMR region

# CHD8 DMR coordinates: chr14:52243603-52243703

CHD8_DMR_region <- subset(dmrs, dmrs$start == "52243603")
plotRegion(merged_BSseq_smooted_filtered, CHD8_DMR_region, extend = 5000, addRegions = CHD8_DMR_region)

# HDGF DMR coordinates: chr3:87908595-87908735

HDGF_DMR_region <- subset(dmrs, dmrs$start == "87908595")
plotRegion(merged_BSseq_smooted_filtered, HDGF_DMR_region, extend = 5000, addRegions = HDGF_DMR_region)

# HR DMR coordinates: chr14:70569594-70570004

HR_DMR_region <- subset(dmrs, dmrs$start == "70569594")
plotRegion(merged_BSseq_smooted_filtered, HR_DMR_region, extend = 5000, addRegions = HR_DMR_region)

# KDM5C DMR coordinates: chrX:152233983-152234421

KDM5C_DMR_region <- subset(dmrs, dmrs$start == "152233983")
plotRegion(merged_BSseq_smooted_filtered, KDM5C_DMR_region, extend = 5000, addRegions = KDM5C_DMR_region)

# MORC4 DMR coordinates: chrX:139870471-139871803

MORC4_DMR_region <- subset(dmrs, dmrs$start == "139870471")
plotRegion(merged_BSseq_smooted_filtered, MORC4_DMR_region, extend = 5000, addRegions = MORC4_DMR_region)

# PHF8 DMR coordinates: chrX:151499108-151499839

PHF8_DMR_region <- subset(dmrs, dmrs$start == "151499108")
plotRegion(merged_BSseq_smooted_filtered, PHF8_DMR_region, extend = 5000, addRegions = PHF8_DMR_region)

# PRDM5 DMR coordinates: chr6:65933001-65933544

PRDM5_DMR_region <- subset(dmrs, dmrs$start == "65933001")
plotRegion(merged_BSseq_smooted_filtered, PRDM5_DMR_region, extend = 5000, addRegions = PRDM5_DMR_region)

# PWWP2A DMR coordinates: chr11:43683117-43683305

PWWP2A_DMR_region <- subset(dmrs, dmrs$start == "43683117")
plotRegion(merged_BSseq_smooted_filtered, PWWP2A_DMR_region, extend = 5000, addRegions = PWWP2A_DMR_region)

# SCML2 DMR coordinates 1: chrX:161161019-161161113
# SCML2 DMR coordinates 2: chrX:161162080-161162192
# SCML2 DMR coordinates 3: chrX:161162598-161163226

SCML2_DMR_region_1 <- subset(dmrs, dmrs$start == "161161019")
SCML2_DMR_region_2 <- subset(dmrs, dmrs$start == "161162080")
SCML2_DMR_region_3 <- subset(dmrs, dmrs$start == "161162598")

plotRegion(merged_BSseq_smooted_filtered, SCML2_DMR_region_1, extend = 5000, addRegions = SCML2_DMR_region_1)

# SND1 DMR coordinates: chr6:28927798-28929771

SND1_DMR_region_1 <- subset(dmrs, dmrs$start == "28927798")
plotRegion(merged_BSseq_smooted_filtered, SND1_DMR_region_1, extend = 5000, addRegions = SND1_DMR_region_1)

# TDRD3 DMR coordinates 1: chr14:87569312-87570163
# TDRD3 DMR coordinates 2: chr14:87912448-87912713

TDRD3_DMR_region_1 <- subset(dmrs, dmrs$start == "87569312")
TDRD3_DMR_region_2 <- subset(dmrs, dmrs$start == "87912448")
plotRegion(merged_BSseq_smooted_filtered, TDRD3_DMR_region_1, extend = 5000, addRegions = TDRD3_DMR_region_1)

# TET2 DMR coordinates: chr3:133545392-133545539

TET2_DMR_region_1 <- subset(dmrs, dmrs$start == "133545392")
plotRegion(merged_BSseq_smooted_filtered, TET2_DMR_region_1, extend = 5000, addRegions = TET2_DMR_region_1)

# TRIM66 DMR coordinates: chr7:109494582-109495053

TRIM66_DMR_region_1 <- subset(dmrs, dmrs$start == "109494582")
plotRegion(merged_BSseq_smooted_filtered, TRIM66_DMR_region_1, extend = 5000, addRegions = TRIM66_DMR_region_1)

# Showing EM genes' DMRs with heatmap start:

# 1. Subset DMRs according to start positions

subset_DMRs <- subset(dmrs, dmrs$start %in% c("52243603", "87908595", "70569594", "152233983", "139870471", "151499108",
                                              "65933001", "43683117", "161161019", "161162080", "161162598",
                                              "28927798", "87569312", "87912448", "133545392", "109494582"))

# 2. Generation GRanges coordinates for DMRs

GRanges <- GRanges(seqnames = c(as.character(subset_DMRs$chr)),
                   ranges = IRanges(start = subset_DMRs$start, end = subset_DMRs$end))

# 3. Obtaining smoothed methylation values of DMRs

methylation_values <- getMeth(merged_BSseq_smooted_filtered, GRanges, type = "smooth", what = "perRegion")

# 4. Manually adding DMR names to methylation matrix

write.csv(methylation_values, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/smoothed_methylation_values_of_DMRs.csv", row.names = F, quote = F)

# 4.1. get methylation value of one DMR region as GRanges object

TET2_DMR_region_3.1 <- GRanges(seqnames = c(as.character(TET2_DMR_region_1$chr)),
                           ranges = IRanges(start = TET2_DMR_region_1$start, end = TET2_DMR_region_1$end))

# 4.2. Get methylation value of this particular DMR and find methylation values in previous metyhlation matrix to label DMRs!

getMeth(merged_BSseq_smooted_filtered, TET2_DMR_region_3.1, type = "smooth", what = "perRegion")

# 5. Uploading modified metylation matrix

modified_methylation_matrix <- read.csv("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/DMRs_and_nearest_genes/EM_genes_with_DMRs/smoothed_methylation_values_of_DMRs.csv", row.names = 1)

# 5. Drawing heatmap based on methylation values

pheatmap::pheatmap(modified_methylation_matrix, color=colorRampPalette(c("navy", "gray90", "red"))(50), cellwidth = 50, cellheight = 30)

# Showing EM genes' DMRs with heatmap end

# --------------------- Kmt2a end ---------------------

# --- plotting gene/promoter/enhancer tracks with plotRegion() function in bsseq R package start ---

# Convert bed file to GRanges object

gene_coordinates <- "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_for_Kmt2a_and_controls/UCSC_mm10_gene_coordinates/mm10.refGene.sorted.bed"
gene_coordinates <- import(gene_coordinates)
gene_coordinates

CHD8_intron_1_coordinates <- subset(gene_coordinates, gene_coordinates$name == "NM_201637")
CHD8_intron_1_coordinates <- GRangesList(intron_1 = CHD8_intron_1_coordinates)

# Yukarıdaki kısımda farklı farklı coordinate'ler ayrı sekilde Granges object olarak atandıktan sonra GRangesList ile
# birlestirilip plot'ta annote edilebilir. Ornegin promoter, TFBS, enhancer coordinate' ler vs..

# Gene track is not working!

Chd8_gene_coordinates <- data.frame("chr" = "chr14", "start" =52198151, "end" =52257760, "gene_ID" ="ENSMUSG00000053754", "exon_number"=38,"strand"="-", "gene_name"="Chd8", "isoforms" = c("ENSMUST00000200169"))

# CHD8 DMR coordinates: chr14:52243603-52243703

CHD8_DMR_region <- subset(dmrs, dmrs$start == "52243603")
plotRegion(merged_BSseq_smooted_filtered, CHD8_DMR_region, extend = 5000, addRegions = CHD8_DMR_region, main = "CHD8", annoTrack = CHD8_intron_1_coordinates, cex.anno = 1.3, lwd = 2.5, mainWithWidth = F, regionCol = "#1b9e77")

# --- plotting gene/promoter/enhancer tracks with plotRegion() function in bsseq R package end ---




# Visualization of methylation data with the R package EnrichedHeatmap

# Rationale: By default every genomic signal (methylation, ChIP-seq, or expression data)
# tries to intersect to every target region (coordinates of genes, promoters, CpG islands etc).

# Enriched heatmap is a special type of heatmap which visualizes the enrichment of genomic signals over 
# specific target regions. It is broadly used to visualize e.g. how histone modifications are enriched at 
# transcription start sites (TSS).

# Loading libraries

library(bsseq)
library(EnrichedHeatmap)

# Get signal data start ------

# Loading methylation data in BSseq format

BSseq <- readRDS("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_all_KOs_versus_all_control_groups/Last_analysis_with_last_pipeline/merged_BSseq_smoothed_and_signal_noise_filtered.rds")

# Obtain first sample as BSseq object

colData(BSseq) # First sample is Ctr6_1
rownames(colData(BSseq))

# Subset first sample as BSseq object

Ctr6_1 <- BSseq[,1]
colData(Ctr6_1) # check whether correct sample is obtained or not

# Generate GRanges object for this sample

Ctr6_1_GRanges <- granges(Ctr6_1)

# Add smoothed methylation values to GRanges object

Ctr6_1_GRanges$methylation <- getMeth(Ctr6_1, type = "smooth")

# Get signal data end ------

# Get target regions start ------

# Loading RefSeq gene coordinates of mm10 mouse genome in bed file format

RefSeq_genes <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Enriched_heatmap_visualization/UCSC_mm10_genome.refGene.gene.coordinates.sorted.bed")

# Keep only gene coordinates, gene names, and strand information

RefSeq_genes <- RefSeq_genes[,c(1:4, 6)]

# Add column names to object

colnames(RefSeq_genes) <- c("chromosome", "start", "end", "RefSeqID", "strand")

# Convert RefSeq gene coordinates in bed file to a GRanges object

RefSeq_genes_GRanges <- GRanges(seqnames = RefSeq_genes$chromosome, ranges = IRanges(start = RefSeq_genes$start,
                                                                                     end = RefSeq_genes$end), strand = RefSeq_genes$strand, RefSeqIDs = RefSeq_genes$RefSeqID)

# We first visualize how methylation is enriched at transcription start sites (TSS).

# We extract TSS of genes (note tss has strand information):

tss = promoters(RefSeq_genes_GRanges, upstream = 0, downstream = 1)

# Get target regions end ------

# Similar as other tools, the task of visualization are separated into two steps:

# 1. get association between genomic signals and targets by normalizing into a matrix.
# 2. visualize the matrix by heatmap.

methylation_matrix = normalizeToMatrix(Ctr6_1_GRanges, tss, value_column = "methylation", mean_mode = "absolute",
                         extend = 5000, w = 50, background = NA)

# Assign colors for methylation values

methylation_color = colorRampPalette(c("navy", "gray90", "red"))(50)

# Draw heatmap showing methylation over transcription start sites

EnrichedHeatmap(methylation_matrix, col = methylation_color,
                name = "methylation", column_title = "methylation near TSS",
                top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

# mm10 mouse genome CpG islands as target region

CpG_islands <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Enriched_heatmap_visualization/UCSC_mm10_mouse_genome_CpG_islands.sorted.bed")
colnames(CpG_islands) <- c("chromosomes", "start", "end", "CpGislands")
  
# Generate CpG island GRanges object

CpG_islands_GRanges <- GRanges(seqnames = CpG_islands$chromosomes, ranges = IRanges(start = CpG_islands$start,
                                                                                     end = CpG_islands$end), CpG_island = CpG_islands$CpGislands)

methylation_matrix_2 = normalizeToMatrix(Ctr6_1_GRanges, CpG_islands_GRanges, value_column = "methylation", mean_mode = "absolute",
                                       extend = 5000, w = 50, background = NA, target_ratio = 0.3)

EnrichedHeatmap(methylation_matrix_2, col = methylation_color,
                name = "methylation", column_title = "methylation near CpG islands",
                top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

# Draw heatmap showing methylation signals over only two target regions:

two_CpG_islands <- CpG_islands_GRanges[1:2,] # Fetch out only two CpG islands

methylation_matrix_3 = normalizeToMatrix(Ctr6_1_GRanges, two_CpG_islands, value_column = "methylation", mean_mode = "absolute",
                                         extend = 5000, w = 50, background = NA)

EnrichedHeatmap(methylation_matrix_3, col = methylation_color,
                row_split = sample(c("CpG island 1", "CpG island 2"), length(two_CpG_islands), replace = TRUE),
                name = "methylation", column_title = "methylation near CpG islands",
                top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    ))))

# how to draw methlation pattern of multiple samples for a target region? ---- start

# Ctr6_1 start ----

Ctr6_1 <- BSseq[,1]

# Generate GRanges object for this sample

Ctr6_1_GRanges <- granges(Ctr6_1)

# Add smoothed methylation values

Ctr6_1_GRanges$methylation <- getMeth(Ctr6_1, type = "smooth")

# Ctr6_1 end ----

# Ctr6_2 start ----

Ctr6_2 <- BSseq[,2]

# Generate GRanges object for this sample

Ctr6_2_GRanges <- granges(Ctr6_2)

# Add smoothed methylation values

Ctr6_2_GRanges$methylation <- getMeth(Ctr6_2, type = "smooth")

# Ctr6_2 end ----

# Ctr1 start ----

Ctr1 <- BSseq[,3]

# Generate GRanges object for this sample

Ctr1_GRanges <- granges(Ctr1)

# Add smoothed methylation values

Ctr1_GRanges$methylation <- getMeth(Ctr1, type = "smooth")

# Ctr1 end ----

# Ctr2 start ----

Ctr2 <- BSseq[,4]

# Generate GRanges object for this sample

Ctr2_GRanges <- granges(Ctr2)

# Add smoothed methylation values

Ctr2_GRanges$methylation <- getMeth(Ctr2, type = "smooth")

# Ctr2 end ----

# Ctr3 start ----

Ctr3 <- BSseq[,5]

# Generate GRanges object for this sample

Ctr3_GRanges <- granges(Ctr3)

# Add smoothed methylation values

Ctr3_GRanges$methylation <- getMeth(Ctr3, type = "smooth")

# Ctr3 end ----

# Kmt2a_1 start ----

Kmt2a_1 <- BSseq[,6]

# Generate GRanges object for this sample

Kmt2a_1_GRanges <- granges(Kmt2a_1)

# Add smoothed methylation values

Kmt2a_1_GRanges$methylation <- getMeth(Kmt2a_1, type = "smooth")

# Kmt2a_1 end ----

# Kmt2a_2 start ----

Kmt2a_2 <- BSseq[,7]

# Generate GRanges object for this sample

Kmt2a_2_GRanges <- granges(Kmt2a_2)

# Add smoothed methylation values

Kmt2a_2_GRanges$methylation <- getMeth(Kmt2a_2, type = "smooth")

# Kmt2a_2 end ----

# Chd1 start ----

Chd1 <- BSseq[,8]

# Generate GRanges object for this sample

Chd1_GRanges <- granges(Chd1)

# Add smoothed methylation values

Chd1_GRanges$methylation <- getMeth(Chd1, type = "smooth")

# Chd1 end ----

# Crebbp start ----

Crebbp <- BSseq[,9]

# Generate GRanges object for this sample

Crebbp_GRanges <- granges(Crebbp)

# Add smoothed methylation values

Crebbp_GRanges$methylation <- getMeth(Crebbp, type = "smooth")

# Crebbp end ----

# Dnmt1 start ----

Dnmt1 <- BSseq[,10]

# Generate GRanges object for this sample

Dnmt1_GRanges <- granges(Dnmt1)

# Add smoothed methylation values

Dnmt1_GRanges$methylation <- getMeth(Dnmt1, type = "smooth")

# Dnmt1 end ----

# Ezh2 start ----

Ezh2 <- BSseq[,11]

# Generate GRanges object for this sample

Ezh2_GRanges <- granges(Ezh2)

# Add smoothed methylation values

Ezh2_GRanges$methylation <- getMeth(Ezh2, type = "smooth")

# Ezh2 end ----

# Hdac6 start ----

Hdac6 <- BSseq[,12]

# Generate GRanges object for this sample

Hdac6_GRanges <- granges(Hdac6)

# Add smoothed methylation values

Hdac6_GRanges$methylation <- getMeth(Hdac6, type = "smooth")

# Hdac6 end ----

# Hdac8 start ----

Hdac8 <- BSseq[,13]

# Generate GRanges object for this sample

Hdac8_GRanges <- granges(Hdac8)

# Add smoothed methylation values

Hdac8_GRanges$methylation <- getMeth(Hdac8, type = "smooth")

# Hdac8 end ----

# Kdm1a start ----

Kdm1a <- BSseq[,14]

# Generate GRanges object for this sample

Kdm1a_GRanges <- granges(Kdm1a)

# Add smoothed methylation values

Kdm1a_GRanges$methylation <- getMeth(Kdm1a, type = "smooth")

# Kdm1a end ----

# Kdm2b start ----

Kdm2b <- BSseq[,15]

# Generate GRanges object for this sample

Kdm2b_GRanges <- granges(Kdm2b)

# Add smoothed methylation values

Kdm2b_GRanges$methylation <- getMeth(Kdm2b, type = "smooth")

# Kdm2b end ----

# Kdm5b start ----

Kdm5b <- BSseq[,16]

# Generate GRanges object for this sample

Kdm5b_GRanges <- granges(Kdm5b)

# Add smoothed methylation values

Kdm5b_GRanges$methylation <- getMeth(Kdm5b, type = "smooth")

# Kdm5b end ----

# Kdm5c start ----

Kdm5c <- BSseq[,17]

# Generate GRanges object for this sample

Kdm5c_GRanges <- granges(Kdm5c)

# Add smoothed methylation values

Kdm5c_GRanges$methylation <- getMeth(Kdm5c, type = "smooth")

# Kdm5c end ----

# Kdm6a start ----

Kdm6a <- BSseq[,18]

# Generate GRanges object for this sample

Kdm6a_GRanges <- granges(Kdm6a)

# Add smoothed methylation values

Kdm6a_GRanges$methylation <- getMeth(Kdm6a, type = "smooth")

# Kdm6a end ----

# Create normalized matrix for each sample

matrix_Ctr6_1 = normalizeToMatrix(Ctr6_1_GRanges, tss, value_column = "methylation", 
                         extend = 5000, mean_mode = "absolute", w = 50, background = NA)

#

matrix_Ctr6_2 = normalizeToMatrix(Ctr6_2_GRanges, tss, value_column = "methylation", 
                                  extend = 5000, mean_mode = "absolute", w = 50, background = NA)

#

matrix_Kmt2a_1 = normalizeToMatrix(Kmt2a_1_GRanges, tss, value_column = "methylation", 
                                  extend = 5000, mean_mode = "absolute", w = 50, background = NA)

#

matrix_Kmt2a_2 = normalizeToMatrix(Kmt2a_2_GRanges, tss, value_column = "methylation", 
                                   extend = 5000, mean_mode = "absolute", w = 50, background = NA)

# Visualize both Ctr6_1 and Ctr6_2 control samples' methylation values over target TSS!
# To combine two samples' methylation values over target sites add "+" character between EnrichedHeatmap()
# functions

EnrichedHeatmap(matrix_Ctr6_1, col = methylation_color,
                name = "methylation", column_title = "Ctr6_1",
                top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    )))) + 
  EnrichedHeatmap(matrix_Ctr6_2, col = methylation_color,
                  name = "methylation", column_title = "Ctr6_2",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      )))) + 
  EnrichedHeatmap(matrix_Kmt2a_1, col = methylation_color,
                  name = "methylation", column_title = "Kmt2a_1",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      )))) + 
  EnrichedHeatmap(matrix_Kmt2a_2, col = methylation_color,
                  name = "methylation", column_title = "Kmt2a_2",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      ))))

# How to draw methlation pattern of multiple samples for a target region? ---- end

# Draw methylation values over vicinity of TSS for neuron-specific genes ----- start

# Neuron-specific genes were downloaded from: "https://doi.org/10.1016/j.celrep.2017.08.086"

Neuron_specific_genes <- read_excel(path = "/home/ko/Downloads/1-s2.0-S2211124717312287-mmc2.xls", sheet = 1, skip = 1)

# mm10 mouse genome RefSeq genes with gene symbols

Refseq_genes_with_gene_symbols <- fread("/home/ko/Downloads/Hans_thought_about_the_most_dramatic_samples_in_terms_of_methylation_difference_from_control_samples/RefSeq_mouse_gene_IDs.csv")

# To get mm10 coordinates of neuron-specific genes from RefSeq

Neuron_specific_genes_with_coordinates <- subset(Refseq_genes_with_gene_symbols, Refseq_genes_with_gene_symbols$GeneSymbols %in% Neuron_specific_genes$`Official Symbol`)

# Generate GRanges object for neuron-specific genes

Neuron_specific_genes_GRanges <- GRanges(seqnames = Neuron_specific_genes_with_coordinates$chromosome, ranges = IRanges(start = Neuron_specific_genes_with_coordinates$start,
                                                                                                                        end = Neuron_specific_genes_with_coordinates$end), Genes = Neuron_specific_genes_with_coordinates$GeneSymbols)
# Extract TSS of neuron-specific genes

Neuron_specific_genes_GRanges_tss = promoters(Neuron_specific_genes_GRanges, upstream = 0, downstream = 1)

# How Dnmt1 knockout methylation values over TSS of neuron specific genes compare to control samples

matrix_Dnmt1 = normalizeToMatrix(Dnmt1_GRanges, Neuron_specific_genes_GRanges_tss, value_column = "methylation", 
                                   extend = 5000, mean_mode = "absolute", w = 50, background = NA)

matrix_Ctr1 = normalizeToMatrix(Ctr1_GRanges, Neuron_specific_genes_GRanges_tss, value_column = "methylation", 
                                 extend = 5000, mean_mode = "absolute", w = 50, background = NA)

matrix_Ctr2 = normalizeToMatrix(Ctr2_GRanges, Neuron_specific_genes_GRanges_tss, value_column = "methylation", 
                                extend = 5000, mean_mode = "absolute", w = 50, background = NA)

matrix_Ctr3 = normalizeToMatrix(Ctr3_GRanges, Neuron_specific_genes_GRanges_tss, value_column = "methylation", 
                                extend = 5000, mean_mode = "absolute", w = 50, background = NA)

# Visualization of global methylation in Dnmt1 and control samples over transcription start sites of
# neuron-specific genes

# Specify the file path and name for the PDF file

pdf_file <- "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Enriched_heatmap_visualization/neuron_specific_genes_TSS_methylation_values.pdf"

# Open the PDF device
pdf(pdf_file)

ht_list <- EnrichedHeatmap(matrix_Dnmt1, col = methylation_color,
                name = "methylation", column_title = "Dnmt1",
                top_annotation = HeatmapAnnotation(
                  enriched = anno_enriched(
                    ylim = c(0, 1),
                    axis_param = list(
                      at = c(0, 0.5, 1),
                      labels = c("0", "0.5", "1"),
                      side = "right",
                      facing = "outside"
                    )))) + 
  EnrichedHeatmap(matrix_Ctr1, col = methylation_color,
                  name = "methylation", column_title = "Ctr1",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      )))) + 
  EnrichedHeatmap(matrix_Ctr2, col = methylation_color,
                  name = "methylation", column_title = "Ctr2",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      )))) + 
  EnrichedHeatmap(matrix_Ctr3, col = methylation_color,
                  name = "methylation", column_title = "Ctr3",
                  top_annotation = HeatmapAnnotation(
                    enriched = anno_enriched(
                      ylim = c(0, 1),
                      axis_param = list(
                        at = c(0, 0.5, 1),
                        labels = c("0", "0.5", "1"),
                        side = "right",
                        facing = "outside"
                      ))))

draw(ht_list, ht_gap = unit(c(8, 8, 8), "mm")) # draw() function put space between heatmaps

# Close the PDF device
dev.off()

# Draw methylation values over vicinity of TSS for neuron-specific genes ----- end

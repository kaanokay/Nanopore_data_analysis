# Visualization of UCSC GREAT enrichment results of DMRs

# --- Start ---

# Input is UCSC GREAT's output!

# GREAT's link: "http://great.stanford.edu/great/public-4.0.4/html/"

# Uploading libraries

library(ggplot2)
library(data.table)

# Uploading UCSC GREAT output

Kmt2a_392_DMR_enrichment_results <- fread("/home/ko/Downloads/GREAT_enrichment_of_DMRs/input_392_DMR_coordinates/outputs/gene_matrices/all_enrichment_results_for_392_DMRs.csv", header = F)
colnames(Kmt2a_392_DMR_enrichment_results) <- c("Ontology", "ID", "Desc", "BinomRank", "BinomP", "BinomBonfP", "BinomFdrQ",
                                                "RegionFoldEnrich", "ExpRegions", "ObsRegions", "GenomeFrac", "SetCov",
                                                "HyperRank", "HyperP", "HyperBonfP", "HyperFdrQ", "GeneFoldEnrich",
                                                "ExpGenes", "ObsGenes", "TotalGenes", "GeneSetCov", "TermCov", "Regions",
                                                "Genes")

# GO Biological Process visualization

# Subset of GO Biological Process term

GO_Biological_Process <- subset(Kmt2a_392_DMR_enrichment_results, Kmt2a_392_DMR_enrichment_results$Ontology == "GO Biological Process")

# Put cutoff for adjusted pvalue and enrichment fold change

GO_Biological_Process <- subset(GO_Biological_Process, GO_Biological_Process$BinomFdrQ < 0.05 & GO_Biological_Process$RegionFoldEnrich > 2)
# subset of enrichments according to adjusted pvalues and fold enrichment coefficient

# Why we used BinomFdrQ and Fold change cutoffs: look into following paper "https://www.nature.com/articles/nature14248)"
# We restricted ourselves to interpretation of results with an enrichment ratio of at least 2,
# and multiple hypothesis testing corrected P values <0.01 (example for GREAT cutoff

# Assignment of fold enrichment as numeric to get top 20 the most enriched terms

GO_Biological_Process$RegionFoldEnrich <- as.numeric(GO_Biological_Process$RegionFoldEnrich)

# Select top 20 enrichment for RegionFoldEnrich column (top 20 enrichment foldchange)

decimal_column <- GO_Biological_Process$RegionFoldEnrich
top_20_enrichment <- which(decimal_column %in% tail(sort(decimal_column), 20))
GO_Biological_Process_top_20_enrichment <- GO_Biological_Process[top_20_enrichment, , drop = FALSE]
print(GO_Biological_Process_top_20_enrichment)

# Selection of Ontology terms, adjusted pvalues and Fold enrichments

GO_Biological_Process_top_20_enrichment <- GO_Biological_Process_top_20_enrichment[,c(3,7,8)]
colnames(GO_Biological_Process_top_20_enrichment) <- c("Term", "FDR", "FoldChange")
GO_Biological_Process_top_20_enrichment$FDR <- as.numeric(GO_Biological_Process_top_20_enrichment$FDR)
GO_Biological_Process_top_20_enrichment$FoldChange <- as.numeric(GO_Biological_Process_top_20_enrichment$FoldChange)

# Sort the data by FoldChange in ascending order
GO_Biological_Process_top_20_enrichment <- GO_Biological_Process_top_20_enrichment[order(GO_Biological_Process_top_20_enrichment$FoldChange),]

# Create a factor variable for the Term column to preserve the desired order
GO_Biological_Process_top_20_enrichment$Term <- factor(GO_Biological_Process_top_20_enrichment$Term, levels = GO_Biological_Process_top_20_enrichment$Term)

# Create a reversed color scale for adjusted p-values

rev_pvalue_scale <- function() {
  trans_new("reverse", function(x) 1 - x, function(x) 1 - x, domain = c(0, 1))
}

# Create the bubble plot
bubble_plot <- ggplot(GO_Biological_Process_top_20_enrichment, aes(x = FoldChange, y = Term, size = FoldChange, color = FDR)) +
  geom_point() +
  scale_size(range = c(3, 10)) +
  scale_color_gradientn(colors = c("blue", "red"), trans = rev_pvalue_scale(), guide = "colorbar", labels = comma_format()) +
  labs(x = "Fold Change", y = NULL) +
  theme_minimal()

# Adjust plot appearance
bubble_plot <- bubble_plot +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# Save the plot as a PDF file
ggsave("GO_Biological_Process.pdf", plot = bubble_plot, device = "pdf", width = 11.7, height = 8.3, units = "in")

# Confirmation message
cat("Bubble plot saved as 'GO_Biological_Process.pdf' in the current directory.\n")

# --- End ---
# --- Fisher test start ---

# Uploading libraries

library(bsseq)
library(data.table)

BSseq.obj <- readRDS("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_all_KOs_versus_all_control_groups/Last_analysis_with_last_pipeline/merged_BSseq_smoothed_filtered.rds")

# Chd1 fisher exact test

Chd1_fisher_test <- fisherTests(BSseq.obj, group1 = "Chd1", group2 = c("Ctr6_1", "Ctr6_2", "Ctr1", "Ctr2", "Ctr3"), mc.cores = 8)

# Obtaining results as a dataframe

Chd1_fisher_results <- as.data.frame(Chd1_fisher_test$results)

# Multiple correction of pvalues with bonferroni test

Chd1_fisher_results$Adjusted_pvalues = p.adjust(Chd1_fisher_results$p.value, method = "bonferroni")

# Obtaining statistically significant DMRs

Chd1_significant_CpGs <- subset(Chd1_fisher_results,
                                Chd1_fisher_results$Adjusted_pvalues < 0.05)

# In Chd1_significant_CpGs, each row represents which row is differentially methylated in BSseq object.

# To get coordinates of statistically significant DMRs from fisher test:

Chd1_significant_CpG_coordinates <- granges(BSseq.obj)[which(Chd1_fisher_results$Adjusted_pvalues < 0.05) ,]

# Adding coordinates to DMR dataframe

chromosome <- as.data.frame(Chd1_significant_CpG_coordinates@seqnames)
coordinates <- as.data.frame(Chd1_significant_CpG_coordinates@ranges)

# Combine significant CpG dataframe with coordinate dataframes

Chd1_significant_CpGs <- cbind(Chd1_significant_CpGs, chromosome)
Chd1_significant_CpGs <- cbind(Chd1_significant_CpGs, coordinates)

colnames(Chd1_significant_CpGs) <- c("p.value", "log2OR", "Adjusted_pvalues", "chromosome", "start", "end", "width")
write.csv(Chd1_significant_CpGs, "/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Fisher_exact_test_results/Chd1_significant_CpGs_fisher_exact_test.csv", quote = F, row.names = F, sep = "\t")

# --- Fisher test end ---
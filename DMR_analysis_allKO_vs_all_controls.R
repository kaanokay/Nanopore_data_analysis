# Loading libraries

library(rtracklayer)
library(bsseq)
library(BiocParallel)
library(parallel)
library(tidyr)
library(tidyselect)
library(RColorBrewer)
library(viridisLite)
library(viridis)
library(scales)
library(colorspace)
library(BRGenomics)
library(data.table)

# Reading modbam2bed outputs

Kmt2a_1_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Kmt2a_experiment/Kmt2a_1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Kmt2a_2_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Kmt2a_experiment/Kmt2a_2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Chd1_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Chd1_experiment/Chd1_CpGs_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Crebbp_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Crebbp_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Dnmt1_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Dnmt1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Ezh2_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Ezh2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Hdac6_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Hdac6_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Hdac8_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Kdm1a_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Kdm1a_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Kdm2b_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Kdm2b_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Kdm5b_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Kdm5b_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Kdm5c_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Kdm5c_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Kdm6a_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Kdm6a_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Ctr6_1_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Kmt2a_experiment/Ctr6_1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Ctr6_2_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Kmt2a_experiment/Ctr6_2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Ctr1_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230420_batch_1_experiment/Ctr1_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")
Ctr2_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/Chd1_experiment/Ctr2_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")
Ctr3_bed_file <- fread("/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/20230307_batch_1_experiment/Ctr3_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")

# Assignment of column names

colnames(Kmt2a_1_bed_file) <- c("chrom","start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")
colnames(Kmt2a_2_bed_file) <- c("chrom","start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")
colnames(Chd1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")
colnames(Crebbp_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Dnmt1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Ezh2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Hdac6_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Hdac8_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Kdm1a_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Kdm2b_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Kdm5b_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Kdm5c_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Kdm6a_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Ctr6_1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                               "freq", "canon", "mod", "filt")
colnames(Ctr6_2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                               "freq", "canon", "mod", "filt")
colnames(Ctr1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Ctr2_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
colnames(Ctr3_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                             "freq", "canon", "mod", "filt")
                             

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kmt2a_1

Kmt2a_1_bed_file$mod_and_canon <- Kmt2a_1_bed_file$mod + Kmt2a_1_bed_file$canon

# BSseq object for Kmt2a_1

Kmt2a_1_BSseq <- BSseq(chr = as.character(Kmt2a_1_bed_file$chrom), pos = Kmt2a_1_bed_file$start +1, 
                       M = matrix(Kmt2a_1_bed_file$mod), Cov = matrix(Kmt2a_1_bed_file$mod_and_canon), sampleNames = c("Kmt2a_1"))
sampleNames(Kmt2a_1_BSseq) <- "Kmt2a_1"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kmt2a_2

Kmt2a_2_bed_file$mod_and_canon <- Kmt2a_2_bed_file$mod + Kmt2a_2_bed_file$canon

# BSseq object for Kmt2a_2

Kmt2a_2_BSseq <- BSseq(chr = as.character(Kmt2a_2_bed_file$chrom), pos = Kmt2a_2_bed_file$start +1, 
                       M = matrix(Kmt2a_2_bed_file$mod), Cov = matrix(Kmt2a_2_bed_file$mod_and_canon), sampleNames = c("Kmt2a_2"))
sampleNames(Kmt2a_2_BSseq) <- "Kmt2a_2"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Chd1

Chd1_bed_file$mod_and_canon <- Chd1_bed_file$mod + Chd1_bed_file$canon

# BSseq object for Chd1

Chd1_BSseq <- BSseq(chr = as.character(Chd1_bed_file$chrom), pos = Chd1_bed_file$start +1, 
                       M = matrix(Chd1_bed_file$mod), Cov = matrix(Chd1_bed_file$mod_and_canon), sampleNames = c("Chd1"))
sampleNames(Chd1_BSseq) <- "Chd1"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Crebbp

Crebbp_bed_file$mod_and_canon <- Crebbp_bed_file$mod + Crebbp_bed_file$canon

# BSseq object for Crebbp

Crebbp_BSseq <- BSseq(chr = as.character(Crebbp_bed_file$chrom), pos = Crebbp_bed_file$start +1, 
                       M = matrix(Crebbp_bed_file$mod), Cov = matrix(Crebbp_bed_file$mod_and_canon), sampleNames = c("Crebbp"))
sampleNames(Crebbp_BSseq) <- "Crebbp"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Dnmt1

Dnmt1_bed_file$mod_and_canon <- Dnmt1_bed_file$mod + Dnmt1_bed_file$canon

# BSseq object for Dnmt1

Dnmt1_BSseq <- BSseq(chr = as.character(Dnmt1_bed_file$chrom), pos = Dnmt1_bed_file$start +1, 
                       M = matrix(Dnmt1_bed_file$mod), Cov = matrix(Dnmt1_bed_file$mod_and_canon), sampleNames = c("Dnmt1"))
sampleNames(Dnmt1_BSseq) <- "Dnmt1"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ezh2

Ezh2_bed_file$mod_and_canon <- Ezh2_bed_file$mod + Ezh2_bed_file$canon

# BSseq object for Ezh2

Ezh2_BSseq <- BSseq(chr = as.character(Ezh2_bed_file$chrom), pos = Ezh2_bed_file$start +1, 
                       M = matrix(Ezh2_bed_file$mod), Cov = matrix(Ezh2_bed_file$mod_and_canon), sampleNames = c("Ezh2"))
sampleNames(Ezh2_BSseq) <- "Ezh2"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Hdac6

Hdac6_bed_file$mod_and_canon <- Hdac6_bed_file$mod + Hdac6_bed_file$canon

# BSseq object for Hdac6

Hdac6_BSseq <- BSseq(chr = as.character(Hdac6_bed_file$chrom), pos = Hdac6_bed_file$start +1, 
                       M = matrix(Hdac6_bed_file$mod), Cov = matrix(Hdac6_bed_file$mod_and_canon), sampleNames = c("Hdac6"))
sampleNames(Hdac6_BSseq) <- "Hdac6"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Hdac8

Hdac8_bed_file$mod_and_canon <- Hdac8_bed_file$mod + Hdac8_bed_file$canon

# BSseq object for Hdac8

Hdac8_BSseq <- BSseq(chr = as.character(Hdac8_bed_file$chrom), pos = Hdac8_bed_file$start +1, 
                       M = matrix(Hdac8_bed_file$mod), Cov = matrix(Hdac8_bed_file$mod_and_canon), sampleNames = c("Hdac8"))
sampleNames(Hdac8_BSseq) <- "Hdac8"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kdm1a

Kdm1a_bed_file$mod_and_canon <- Kdm1a_bed_file$mod + Kdm1a_bed_file$canon

# BSseq object for Kdm1a

Kdm1a_BSseq <- BSseq(chr = as.character(Kdm1a_bed_file$chrom), pos = Kdm1a_bed_file$start +1, 
                       M = matrix(Kdm1a_bed_file$mod), Cov = matrix(Kdm1a_bed_file$mod_and_canon), sampleNames = c("Kdm1a"))
sampleNames(Kdm1a_BSseq) <- "Kdm1a"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kdm2b

Kdm2b_bed_file$mod_and_canon <- Kdm2b_bed_file$mod + Kdm2b_bed_file$canon

# BSseq object for Kdm2b

Kdm2b_BSseq <- BSseq(chr = as.character(Kdm2b_bed_file$chrom), pos = Kdm2b_bed_file$start +1, 
                       M = matrix(Kdm2b_bed_file$mod), Cov = matrix(Kdm2b_bed_file$mod_and_canon), sampleNames = c("Kdm2b"))
sampleNames(Kdm2b_BSseq) <- "Kdm2b"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kdm5b

Kdm5b_bed_file$mod_and_canon <- Kdm5b_bed_file$mod + Kdm5b_bed_file$canon

# BSseq object for Kdm5b

Kdm5b_BSseq <- BSseq(chr = as.character(Kdm5b_bed_file$chrom), pos = Kdm5b_bed_file$start +1, 
                       M = matrix(Kdm5b_bed_file$mod), Cov = matrix(Kdm5b_bed_file$mod_and_canon), sampleNames = c("Kdm5b"))
sampleNames(Kdm5b_BSseq) <- "Kdm5b"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kdm5c

Kdm5c_bed_file$mod_and_canon <- Kdm5c_bed_file$mod + Kdm5c_bed_file$canon

# BSseq object for Kdm5c

Kdm5c_BSseq <- BSseq(chr = as.character(Kdm5c_bed_file$chrom), pos = Kdm5c_bed_file$start +1, 
                       M = matrix(Kdm5c_bed_file$mod), Cov = matrix(Kdm5c_bed_file$mod_and_canon), sampleNames = c("Kdm5c"))
sampleNames(Kdm5c_BSseq) <- "Kdm5c"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Kdm6a

Kdm6a_bed_file$mod_and_canon <- Kdm6a_bed_file$mod + Kdm6a_bed_file$canon

# BSseq object for Kdm6a

Kdm6a_BSseq <- BSseq(chr = as.character(Kdm6a_bed_file$chrom), pos = Kdm6a_bed_file$start +1, 
                       M = matrix(Kdm6a_bed_file$mod), Cov = matrix(Kdm6a_bed_file$mod_and_canon), sampleNames = c("Kdm6a"))
sampleNames(Kdm6a_BSseq) <- "Kdm6a"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ctr6_1

Ctr6_1_bed_file$mod_and_canon <- Ctr6_1_bed_file$mod + Ctr6_1_bed_file$canon

# BSseq object for Ctr6_1

Ctr6_1_BSseq <- BSseq(chr = as.character(Ctr6_1_bed_file$chrom), pos = Ctr6_1_bed_file$start +1, 
                       M = matrix(Ctr6_1_bed_file$mod), Cov = matrix(Ctr6_1_bed_file$mod_and_canon), sampleNames = c("Ctr6_1"))
sampleNames(Ctr6_1_BSseq) <- "Ctr6_1"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ctr6_2

Ctr6_2_bed_file$mod_and_canon <- Ctr6_2_bed_file$mod + Ctr6_2_bed_file$canon

# BSseq object for Ctr6_2

Ctr6_2_BSseq <- BSseq(chr = as.character(Ctr6_2_bed_file$chrom), pos = Ctr6_2_bed_file$start +1, 
                       M = matrix(Ctr6_2_bed_file$mod), Cov = matrix(Ctr6_2_bed_file$mod_and_canon), sampleNames = c("Ctr6_2"))
sampleNames(Ctr6_2_BSseq) <- "Ctr6_2"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ctr1

Ctr1_bed_file$mod_and_canon <- Ctr1_bed_file$mod + Ctr1_bed_file$canon

# BSseq object for Ctr1

Ctr1_BSseq <- BSseq(chr = as.character(Ctr1_bed_file$chrom), pos = Ctr1_bed_file$start +1, 
                       M = matrix(Ctr1_bed_file$mod), Cov = matrix(Ctr1_bed_file$mod_and_canon), sampleNames = c("Ctr1"))
sampleNames(Ctr1_BSseq) <- "Ctr1"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ctr2

Ctr2_bed_file$mod_and_canon <- Ctr2_bed_file$mod + Ctr2_bed_file$canon

# BSseq object for Ctr2

Ctr2_BSseq <- BSseq(chr = as.character(Ctr2_bed_file$chrom), pos = Ctr2_bed_file$start +1, 
                       M = matrix(Ctr2_bed_file$mod), Cov = matrix(Ctr2_bed_file$mod_and_canon), sampleNames = c("Ctr2"))
sampleNames(Ctr2_BSseq) <- "Ctr2"

# Coverage calculation for each CpG by sum canonical and modified CpGs for Ctr3

Ctr3_bed_file$mod_and_canon <- Ctr3_bed_file$mod + Ctr3_bed_file$canon

# BSseq object for Ctr3

Ctr3_BSseq <- BSseq(chr = as.character(Ctr3_bed_file$chrom), pos = Ctr3_bed_file$start +1, 
                       M = matrix(Ctr3_bed_file$mod), Cov = matrix(Ctr3_bed_file$mod_and_canon), sampleNames = c("Ctr3"))
sampleNames(Ctr3_BSseq) <- "Ctr3"

# Combining BSseq objects

merged_BSseq <- combine(Ctr6_1_BSseq, Ctr6_2_BSseq, Ctr1_BSseq, Ctr2_BSseq, Ctr3_BSseq, Kmt2a_1_BSseq, Kmt2a_2_BSseq,
                        Chd1_BSseq, Crebbp_BSseq, Dnmt1_BSseq, Ezh2_BSseq, Hdac6_BSseq, Hdac8_BSseq, Kdm1a_BSseq,
                        Kdm2b_BSseq, Kdm5b_BSseq, Kdm5c_BSseq, Kdm6a_BSseq)

# Collapsedata

merged_BSseq <- collapseBSseq(merged_BSseq, group = c("Control1", "Control2", "Control3", "Control4", "Control5",
                                                      "Kmt2a_1", "Kmt2a_2", "Chd1", "Crebbp", "Dnmt1", "Ezh2",
                                                      "Hdac6", "Hdac8", "Kdm1a", "Kdm2b", "Kdm5b" ,"Kdm5c", "Kdm6a"))

# --- Smooting data start ---

merged_BSseq_smoothed <- BSmooth(
  BSseq = merged_BSseq, 
  BPPARAM = MulticoreParam(workers = 20), 
  verbose = TRUE)

# --- Computing t-statistics start ---

# Creating metadata for all samples

metadata_2 <- data.frame("Samples" = c("Control", "Control", "Control","Control","Control",
                                       "AllKOs", "AllKOs", "AllKOs", "AllKOs", "AllKOs", "AllKOs", "AllKOs",
                                       "AllKOs", "AllKOs", "AllKOs", "AllKOs", "AllKOs", "AllKOs"))

pData(merged_BSseq_smoothed)$metadata_2 <- metadata_2
pData(merged_BSseq_smoothed)

# remove CpGs with little or no coverage
# we will only keep CpGs where at least 2 knockout samples and at least 2 Control samples have at least 10x in coverage

cov <- getCoverage(merged_BSseq_smoothed)
keepLoci.ex <- which(rowSums(cov[, merged_BSseq_smoothed$metadata_2 == "AllKOs"] >= 10) >= 2 &
                       rowSums(cov[, merged_BSseq_smoothed$metadata_2 == "Control"] >= 10) >= 2)

length(keepLoci.ex) 

merged_BSseq_smoothed_filtered <- merged_BSseq_smoothed[keepLoci.ex,]

# Save BSseq object

saveRDS(merged_BSseq_smoothed_filtered, file = "/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/DMR_analysis/merged_BSseq_smoothed_filtered.rds")

# Computing t-statistics

merged_BSseq_smoothed_filtered.tstat <- BSmooth.tstat(merged_BSseq_smoothed_filtered, 
                                                     group1 = c("Kmt2a_1", "Kmt2a_2", "Chd1", "Crebbp", "Dnmt1", "Ezh2",
                                                                "Hdac6", "Hdac8", "Kdm1a", "Kdm2b", "Kdm5b", "Kdm5c", "Kdm6a"),
                                                     group2 = c("Ctr6_1", "Ctr6_2", "Ctr1", "Ctr2", "Ctr3"),
                                                     estimate.var = "group2",
                                                     local.correct = TRUE,
                                                     verbose = TRUE, k = 21)

# Finding DMRs

dmrs0 <- dmrFinder(merged_BSseq_smoothed_filtered.tstat, cutoff = c(-4.6, 4.6))
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
write.table(dmrs, "/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/DMR_analysis/allKO_vs_all_controls_DMRs.txt", quote = F, row.names = F, sep =  "\t")

# Color control and knockout samples for plotting reason

pData <- pData(merged_BSseq_smoothed_filtered)
pData$col <- c(rep("#d95f02", each=5), rep("#7570b3", each=13)) # #d95f02 is control, #7570b3 is knockout
pData(merged_BSseq_smoothed_filtered) <- pData

# Plot all DMRs as pdf file

pdf(file = "/hpcdata/Mimir/kao25/Juan_Nanopore_data/yeni_analiz/DMR_analysis/allKOs_vs_all_controls_DMRs.pdf", width = 10, height = 5)
plotManyRegions(merged_BSseq_smoothed_filtered, dmrs[], extend = 5000, addRegions = dmrs)
dev.off()

# End          
          

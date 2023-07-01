# Signal noise detection start ---

# Uploading libraries

library(data.table)

# Uploading modbam2bed output bed file

Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed")

colnames(Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed) <- c("chrom","start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                                                      "freq", "canon", "mod", "filt")

# Detection of range of coverage in bed file (min and max coverage)

range(Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$coverage)

# regions with extreme covrage because of signal-noise

chr2_signal_noise_region_1 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 98662555, end = 98667332))

chr2_signal_noise_region_2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 181917000, end = 181932100))

chr9_signal_noise_region <- GRanges(seqnames = "chr9", ranges = IRanges(start = 3000003, end = 3038419))

chr11_signal_noise_region <- GRanges(seqnames = "chr11", ranges = IRanges(start = 3123000, end = 3201000))

chr17_signal_noise_region <- GRanges(seqnames = "chr17", ranges = IRanges(start = 39843000, end = 39849000))

chr18_signal_noise_region <- GRanges(seqnames = "chr18", ranges = IRanges(start = 3004000, end = 3007000))

# to get all genome regions by creating GRanges object from bed file

Hdac8_MQ_30_and_400bp_read_length_granges <- GRanges(seqnames = Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$chrom, ranges = IRanges(start = Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$start, end = Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$end))

# Find number of row of signal noise regions by overlapping GRanges object of CpG-level data and signal-noise region data

chr2_signal_noise_region_1_rows <- findOverlaps(Hdac8_MQ_30_and_400bp_read_length_granges, chr2_signal_noise_region_1)
str(chr2_signal_noise_region_1_rows@from) # This give us order of rows of signal noise regions
# Let's say it gave us 11141036 which is one of the signal noise CpG coordinate
# Lets inspect it on data
# There are 73 CpGs overlapping with this signal noise location!

Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed[11141036,] # This position is chr2:98662793-98662795
# and coverage is 14278! This is extremely high coverage.

# To count how many CpG overlapping with signal noise region, countOverlaps() function can be used:

table(countOverlaps(Hdac8_MQ_30_and_400bp_read_length_granges, chr2_signal_noise_region_1))

# This function returns 0 and 1, 0 means there is no overlap, 1 means there is an overlap.

# Detection of CpGs overlapping with signal noise region in chr9

chr9_signal_noise_region_rows <- findOverlaps(Hdac8_MQ_30_and_400bp_read_length_granges, chr9_signal_noise_region)
str(chr9_signal_noise_region_rows@from) # There 777 CpGs overlapping with signal noise region in chr9

# Remove CpGs overlapping with signal noise regions (In total, 850 CpGs should be removed)

Hdac8_filtered <- Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed[!(rownames(Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed) %in% chr2_signal_noise_region_1_rows@from) &
                                                                           !(rownames(Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed) %in% chr9_signal_noise_region_rows@from), ]

# Check number of CpGs after and before filtering
nrow(Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed) # 21,267,293 CpGs
nrow(Hdac8_filtered) # 21,266,443 CpGs

# 850 CpGs have successfully been filtered out!

# Check mean coverage of signal noise regions

# Define interval of target signal noise region (chr2)

target_chromosome <- "chr2"
target_start <- 98662555
target_end <- 98667332

Hdac8_chr2_signal_noise_regions <- Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed[Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$chrom == target_chromosome & Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$start >= target_start & Hdac8_MQ_30_and_400bp_read_length_filtered.cpg.acc.bed$end <= target_end, ]
mean(Hdac8_chr2_signal_noise_regions$coverage) # mean coverage is 6018!

# Check whether there are any signal noise DMRs in fisher exact test results, if there are some,
# check their pvalue and log2 values by generating to see there are true positive DMRs between knockouts
# and control samples Start ------

# Uploading fisher's exact test results

Dnmt1 <- read.csv("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Fisher_exact_test_results/Dnmt1_Fisher_test_results.csv")

# Specify the signal noise region interval to subset
target_chromosome <- "chr9"
target_start <- 3000003
target_end <- 3038419

# Subset the data based on the signal noise region interval in chr9
Dnmt1_chr9_signal_noise <- Dnmt1[Dnmt1$seqnames == target_chromosome & Dnmt1$start >= target_start & Dnmt1$end <= target_end, ]

# Draw histogram plot of pvalues and log2 fold changes

hist(Dnmt1_chr9_signal_noise$p.value)
hist(Dnmt1_chr9_signal_noise$log2OR)

# Obtaining non signal noise regions in chr9

Dnmt1_chr9_non_signal_noise <- subset(Dnmt1, !Dnmt1$start %in% Dnmt1_chr9_signal_noise$start)

# Draw histogram plot of pvalues and log2 fold changes
hist(Dnmt1_chr9_non_signal_noise$p.value)
hist(Dnmt1_chr9_non_signal_noise$log2OR)

# aim of histogram plot of pvalues is to be sure whether we have true positive DMRs between knockout samples
# and control samples.

# Check whether there are any signal noise DMRs in fisher exact test results, if there are some,
# check their pvalue and log2 values by generating to see there are true positive DMRs between knockouts
# and control samples End ------

# Removal of signal noise regions from bed files:
# ----------------------------
target_chromosome <- "chr2"
target_start <- 98662555
target_end <- 98667332
# ----------------------------
target_chromosome_1 <- "chr2"
target_start_1 <- 181917000
target_end_1 <- 181932100
# ----------------------------
target_chromosome_2 <- "chr9"
target_start_2 <- 3000003
target_end_2 <- 3038419
# ----------------------------
target_chromosome_3 <- "chr11"
target_start_3 <- 3123000
target_end_3 <- 3201000
# ----------------------------
target_chromosome_4 <- "chr17"
target_start_4 <- 39843000
target_end_4 <- 39849000
# ----------------------------
target_chromosome_5 <- "chr18"
target_start_5 <- 3004000
target_end_5 <- 3007000
# ----------------------------

Kmt2a_2_bed_file <- Kmt2a_2_bed_file[!(Kmt2a_2_bed_file$chrom == target_chromosome &
                                         Kmt2a_2_bed_file$start >= target_start &
                                         Kmt2a_2_bed_file$end <= target_end |
                                         Kmt2a_2_bed_file$chrom == target_chromosome_1 &
                                         Kmt2a_2_bed_file$start >= target_start_1 &
                                         Kmt2a_2_bed_file$end <= target_end_1 |
                                         Kmt2a_2_bed_file$chrom == target_chromosome_2 &
                                         Kmt2a_2_bed_file$start >= target_start_2 &
                                         Kmt2a_2_bed_file$end <= target_end_2 |
                                         Kmt2a_2_bed_file$chrom == target_chromosome_3 &
                                         Kmt2a_2_bed_file$start >= target_start_3 &
                                         Kmt2a_2_bed_file$end <= target_end_3 |
                                         Kmt2a_2_bed_file$chrom == target_chromosome_4 &
                                         Kmt2a_2_bed_file$start >= target_start_4 &
                                         Kmt2a_2_bed_file$end <= target_end_4 |
                                         Kmt2a_2_bed_file$chrom == target_chromosome_5 &
                                         Kmt2a_2_bed_file$start >= target_start_5 &
                                         Kmt2a_2_bed_file$end <= target_end_5), ]


# Signal noise detection end ---

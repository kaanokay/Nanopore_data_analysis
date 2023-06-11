# --- Start ---

# Find the most dramatic samples in terms of average methylation difference from control samples.

# Uploading libraries

library(data.table)

# Uploading BSseq object for all KOs vs all controls

BSseq.obj <- readRDS("/home/ko/Documents/ONT_data/all_KOs_vs_all_controls_results/merged_BSseq_smoothed_filtered.rds")

# See order of samples

colnames(BSseq.obj)

# Uploading 24 DMRs between all KOs vs all controls

DMRs <- fread("/home/ko/Documents/ONT_data/all_KOs_vs_all_controls_results/allKOs_vs_all_controls_DMRs.txt")

# Generating GRanges object with coordinates of 24 DMRs

all_DMRs <- GRanges(seqnames = c(as.character(DMRs$chr)),
                   ranges = IRanges(start = DMRs$start, end = DMRs$end))

# Obtaining smoothed methylation values of all DMRs

All_DMRs_methylation_values <- getMeth(BSseq.obj, all_DMRs, type = "smooth", what = "perRegion")

# Average methylation value of first DMR for control samples minus methylation value of first DMR for Kmt2a_1 sample
# Here we can get how Kmt2a_1 is away from control samples in absolute values of average methylation!
# First five columns are control samples, sixth column is Kmt2a_1
# Command below is performed for all KO samples

abs((rowSums(All_DMRs_methylation_values[1,1:5])/5) - (All_DMRs_methylation_values[1,6]))

# After obtaining average methylation difference between control and each KO,
# We have those differences for each DMR as a matrix where each row represents each KO,
# whereas each column represents each DMR

matrixx <- read.csv("/home/ko/Documents/ONT_data/all_KOs_vs_all_controls_results/Smoothed_methylation_values_of_24_shared_DMRs.csv", row.names = 1)

# matrix looks like structure below:

#         First_DMR   Second_DMR
# Kmt2a   0.01        0.5
# Chd1    0.3         0.15

# Apply the reorder_matrix function to obtain the reordered matrix
# Here according to average absolute difference of methylation of KOs from controls
# ordered from descending order to get the most dramatic sample from control in terms of methylation level.

reordered_matrix <- apply(matrixx, 2, function(x) rownames(matrixx)[order(-x)])

# Create an example matrix
matrix_data <- read.csv("/home/ko/Documents/ONT_data/all_KOs_vs_all_controls_results/Order_of_each_KO_relative_methylation_values_to_control_samples_for_each_DMR.csv")

# matrix_data looks like that:

#   First_DMR     Second_DMR
#   Kdm5c         Hdac6
#   Kdm1a         Dnmt1

# matrix_data has columns of each DMR and order of KOs for each DMR
# For First_DMR Kdm5c is the most dramatic samples, where in Second_DMR
# the most dramatic one is Hdac6.

# Find the most repeated KO in each row
# Aim of this step is find which KO are repeatet at most in each line
# to find the most dramatic KO across all DMRs

most_repeated <- apply(matrix_data, 1, function(row) {
  counts <- table(row)
  max_count <- max(counts)
  most_repeated <- names(counts)[counts == max_count]
  return(most_repeated)
})

# Print the most repeated KO for each row
print(most_repeated)

# To find the most dramatic sample differ from control
# calculate average methylation level difference of all DMRs in KOs from control samples at once
# and order KOs this differences

abs(sum(colSums(All_DMRs_methylation_values[,1:5]))/5) - (sum(All_DMRs_methylation_values[,6]))

# --- End ---
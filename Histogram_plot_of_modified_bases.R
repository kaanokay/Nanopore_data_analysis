# ---- Start

# Loading libraries

library(data.table)

# Uploading modbam2bed files

Kdm5c_bed_file <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Kdm5c_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")

Ctr1_bed_file <- fread("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/Ctr1_MQ_30_and_400bp_filtered_CpGs.cpg.acc.bed")

colnames(Kdm5c_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")

colnames(Ctr1_bed_file) <- c("chrom",	"start", "end", "name", "score", "strand", "tstart", "tend","color", "coverage",
                                "freq", "canon", "mod", "filt")
                                
# Draw histogram plot of modified cytosines

plot.new()

# Plot the first histogram with transparent bars

hist(Kdm5c_bed_file$freq, col = rgb(0, 0, 1, alpha = 0.5), xlab = "% of modified cytosine", ylab = "Frequency", main = "% frequency of modified cytosines")

# Overlay the second histogram with transparent bars on top

hist(Ctr1_bed_file$freq, col = rgb(1, 0.5, 0, alpha = 0.5), add = TRUE)

# Get the midpoints of the plot window
xmid <- mean(par()$usr[1:2])
ymid <- mean(par()$usr[3:4])

# Add a legend in the middle of the plot

legend(xmid, ymid, legend = c("Kdm5c", "Ctr1"), fill = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0.5, 0, alpha = 0.5)), 
       bg = "white", box.lwd = 0.5, cex = 0.8)

# Add a title

title(main = "% frequency of modified cytosines")

# ---- End

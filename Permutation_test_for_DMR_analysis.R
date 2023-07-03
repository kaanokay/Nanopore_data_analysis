# Permutation test to get pvalues for DMRs:

# Start ------------

# Loading libraries ------------

library(bsseq)
library(data.table)
library(limma)
library(parallel)
library(permute)
library(here)

# Uploading BSseq object for all KO vs all control design

BSseq.obj <- readRDS("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_all_KOs_versus_all_control_groups/Last_analysis_with_last_pipeline/merged_BSseq_smoothed_filtered.rds")

# Creating design matrix ------------

pData(BSseq.obj)$metadata_2[,1] # This returns annotation of samples, that is, Control and KOs!
# metada_2 annotation is a data frame but we need to have character strings instead of data frame

BSseq.obj$group <- pData(BSseq.obj)$metadata_2[,1] # create "group" column in BSseq object describing samples
# as control and knockouts

# Create a design matrix to assing control samples "1", where knockout samples 0.
# Design matrix is kind of descriptive data for linear model.

design <- model.matrix(~ group, colData(BSseq.obj))
colnames(design) <- gsub("group", "", colnames(design))

# Create a contrast matrix

contrasts <- diag(rep(1, ncol(design)))
rownames(contrasts) <- colnames(design)

# create fstat.pipeline() function ------------

fstat.pipeline <- function(BSseq, design, contrasts, cutoff, fac, nperm = 1000,
                           coef = NULL, maxGap.sd = 10 ^ 8, maxGap.dmr = 300,
                           type = "dmrs", mc.cores = 1) {
  type <- match.arg(type, c("dmrs", "blocks"))
  bstat <- BSmooth.fstat(BSseq = BSseq, design = design,
                         contrasts = contrasts)
  bstat <- smoothSds(bstat)
  bstat <- computeStat(bstat, coef = coef)
  dmrs <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
  if (is.null(dmrs)) {
    stop("No DMRs identified. Consider reducing the 'cutoff' from (",
         paste0(cutoff, collapse = ", "), ")")
  }
  idxMatrix <- permuteAll(nperm, design)
  nullDist <- getNullDistribution_BSmooth.fstat(BSseq = BSseq,
                                                idxMatrix = idxMatrix,
                                                design = design,
                                                contrasts = contrasts,
                                                coef = coef,
                                                cutoff = cutoff,
                                                maxGap.sd = maxGap.sd,
                                                maxGap.dmr = maxGap.dmr,
                                                mc.cores = mc.cores)
  fwer <- getFWER.fstat(null = c(list(dmrs), nullDist), type = type)
  dmrs$fwer <- fwer
  meth <- getMeth(BSseq, dmrs, what = "perRegion")
  meth <- t(apply(meth, 1, function(xx) tapply(xx, fac, mean)))
  dmrs <- cbind(dmrs, meth)
  dmrs$maxDiff <- rowMaxs(meth) - rowMins(meth)
  list(bstat = bstat, dmrs = dmrs, idxMatrix = idxMatrix, nullDist = nullDist)
}

# Create BSmooth.fstat() function ------------

BSmooth.fstat <- function(BSseq, design, contrasts, verbose = TRUE){
  stopifnot(is(BSseq, "BSseq"))
  stopifnot(hasBeenSmoothed(BSseq))
  
  ## if(any(rowSums(getCoverage(BSseq)[, unlist(groups)]) == 0))
  ##     warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
  
  if(verbose) cat("[BSmooth.fstat] fitting linear models ... ")
  ptime1 <- proc.time()
  allPs <- getMeth(BSseq, type = "smooth", what = "perBase",
                   confint = FALSE)
  allPs <- as.array(allPs)
  fit <- lmFit(allPs, design)
  fitC <- contrasts.fit(fit, contrasts)
  ## Need
  ##   fitC$coefficients, fitC$stdev.unscaled, fitC$sigma, fitC$cov.coefficients
  ## actuall just need
  ##   tstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
  ##   rawSds <- fitC$sigma
  ##   cor.coefficients <- cov2cor(fitC$cov.coefficients)
  rawSds <- as.matrix(fitC$sigma)
  cor.coefficients <- cov2cor(fitC$cov.coefficients)
  rawTstats <- fitC$coefficients / fitC$stdev.unscaled / fitC$sigma
  names(dimnames(rawTstats)) <- NULL
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if(verbose) cat(sprintf("done in %.1f sec\n", stime))
  
  parameters <- c(BSseq@parameters,
                  list(design = design, contrasts = contrasts))
  stats <- list(rawSds = rawSds,
                cor.coefficients = cor.coefficients,
                rawTstats = rawTstats)
  out <- BSseqStat(gr = granges(BSseq),
                   stats = stats,
                   parameters = parameters)
  out
}

# Create smoothSds() function ------------

smoothSds <- function(BSseqStat, k = 101, qSd = 0.75, mc.cores = 1,
                      maxGap = 10^8, verbose = TRUE) {
  smoothSd <- function(Sds, k, qSd) {
    k0 <- floor(k/2)
    if(all(is.na(Sds))) return(Sds)
    thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
    addSD <- rep(median(Sds, na.rm = TRUE), k0)
    sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
    sSds
  }
  if(is.null(maxGap))
    maxGap <- BSseqStat@parameters[["maxGap"]]
  if(is.null(maxGap))
    stop("need to set argument 'maxGap'")
  if(verbose) cat("[smoothSds] preprocessing ... ")
  ptime1 <- proc.time()
  clusterIdx <- makeClusters(granges(BSseqStat), maxGap = maxGap)
  ptime2 <- proc.time()
  stime <- (ptime2 - ptime1)[3]
  if(verbose) cat(sprintf("done in %.1f sec\n", stime))
  smoothSds <- do.call("c",
                       mclapply(clusterIdx, function(idx) {
                         rawSds <- getStats(BSseqStat,
                                            what = "rawSds")[idx, ]
                         rawSds <- as.array(rawSds)
                         smoothSd(rawSds, k = k, qSd = qSd)
                       }, mc.cores = mc.cores))
  smoothSds <- as.matrix(smoothSds)
  if("smoothSds" %in% names(getStats(BSseqStat)))
    BSseqStat@stats[["smoothSds"]] <- smoothSds
  else
    BSseqStat@stats <- c(getStats(BSseqStat), list(smoothSds = smoothSds))
  BSseqStat
}

# Create makeClusters() function ------------

makeClusters <- function(hasGRanges, maxGap = 10^8) {
  chrOrder <- as.character(runValue(seqnames(hasGRanges)))
  if(anyDuplicated(chrOrder))
    stop("argument 'hasGRanges' is not properly order")
  grBase <- granges(hasGRanges)
  clusters <- reduce(resize(grBase, width = 2*maxGap + 1, fix = "center"))
  start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
  clusters.sp <- split(clusters, seqnames(clusters))
  stopifnot(all(sapply(clusters.sp, function(cluster.gr) {
    if(length(cluster.gr) <= 1) return(TRUE)
    all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
  }))) # are the clusters ordered within the chromosome? This is probably guranteed
  clusters <- Reduce(c, clusters.sp[chrOrder])
  stopifnot(all(chrOrder == runValue(seqnames(clusters))))
  ov <- findOverlaps(grBase, clusters)
  clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
  names(clusterIdx) <- NULL
  clusterIdx
}

# Create computeStat() function ------------

computeStat <- function(BSseqStat, coef = NULL) {
  stopifnot(is(BSseqStat, "BSseqStat"))
  if(is.null(coef)) {
    coef <- 1:ncol(getStats(BSseqStat, what = "rawTstats"))
  }
  raw_tstats <- getStats(BSseqStat, what = "rawTstats")[, coef, drop = FALSE]
  scaled_sds <- getStats(BSseqStat, what = "rawSds") /
    getStats(BSseqStat, what = "smoothSds")
  scaled_sds_matrix <- matrix(
    rep(scaled_sds, ncol(raw_tstats)),
    ncol = ncol(raw_tstats))
  tstats <- raw_tstats * scaled_sds_matrix
  if(length(coef) > 1) {
    cor.coefficients <- getStats(BSseqStat,
                                 what = "cor.coefficients")[coef,coef]
    # NOTE: classifyTestsF() calls as.matrix(tstats) and so realises this
    #       array
    stat <- as.matrix(
      classifyTestsF(tstats, cor.coefficients, fstat.only = TRUE))
    stat.type <- "fstat"
  } else {
    stat <- tstats
    stat.type <- "tstat"
  }
  if("stat" %in% names(getStats(BSseqStat))) {
    BSseqStat@stats[["stat"]] <- stat
    BSseqStat@stats[["stat.type"]] <- stat.type
  } else {
    BSseqStat@stats <- c(getStats(BSseqStat),
                         list(stat = stat, stat.type = stat.type))
  }
  BSseqStat
}

# Do permutation test ------------

set.seed(20190206)

fstat_pipeline <- fstat.pipeline(
  BSseq = BSseq.obj,
  design = design,
  contrasts = contrasts,
  # NOTE: This uses F-stat, not t-stat, so square the cutoff.
  cutoff = 4.6 ^ 2,
  fac = BSseq.obj$group,
  nperm = 1000,
  coef = NULL,
  type = "dmrs",
  mc.cores = 64)

# Error in fstat.pipeline(BSseq = BSseq.obj, design = design, contrasts = contrasts,  : 
# No DMRs identified. Consider reducing the 'cutoff' from (21.16)
# Kasper said that: "This error tells us that there are no putative DMRs identified. Yet, you find putative DMRs when you run the tstat function by yourself, so something is wrong. 
# This is a different result from saying that none of the 24 putative DMRs pass multiple testing correction"

# Save results as rds file ------------

saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("Downloads", "save", "allKO_vs_allcontrols_DMR_results.rds"))

# "Downloads" and "save" are directories where results are saved.

# Or save like this to get final DMRs:

write.csv(subset(fstat_pipeline$dmrs, n >= 3 & abs(maxDiff) >= 0.1), "/home/ko/Downloads/Hans_thought_about_the_most_dramatic_samples_in_terms_of_methylation_difference_from_control_samples/Kaan/permutation_test_DMRs.csv", quote = F, row.names = F)
                  
# What is difference between t-stat test and f-stat test in bsseq pipeline:
# Kasper said that "One difference is the estimate.var="group2" That option doesn't really exist for the F-stat approach".

# Ask something about permutation test to ChatGPT:

# Question : what is design matrix and contrasts matrix in statistics, explain me diag function in base R package and model.matrix function in stats R package

# Answer:
# In statistics, the design matrix and contrasts matrix are two important concepts used in regression analysis and statistical modeling.

# Design Matrix: The design matrix, often denoted as X, is a matrix that represents the predictor variables in a regression model. Each row of the design matrix corresponds to an observation, and each column represents a different predictor variable or factor level. In a simple linear regression model with one predictor variable, the design matrix would have two columns: one for the intercept term (usually a column of ones) and one for the predictor variable. In more complex models, such as multiple linear regression or generalized linear models, the design matrix will have additional columns for each predictor or factor level.

# Contrasts Matrix: The contrasts matrix is a matrix used to specify the coding scheme for categorical predictor variables or factor levels in a regression model. It determines how the levels of a categorical variable are represented numerically in the design matrix. Different coding schemes can be used to represent categorical variables, such as dummy coding, treatment coding, or contrast coding. The contrasts matrix defines the numerical values assigned to each level of the categorical variable.

# Now let's move on to the functions in R:

# The diag() function in the base R package is used to create or extract a diagonal matrix or diagonal elements from a given matrix. It has different applications depending on how it is used. Here are some examples:

# To create a diagonal matrix from a vector, you can pass the vector as an argument to the diag() function. For example, diag(c(1, 2, 3)) will create a diagonal matrix with the elements 1, 2, and 3 on its diagonal and zeros elsewhere.
# To extract the diagonal elements from a matrix, you can pass the matrix as an argument to the diag() function without specifying any additional parameters. For example, diag(matrix(1:9, nrow = 3)) will return a vector with the diagonal elements 1, 5, and 9.
# The model.matrix() function is part of the stats package in R and is used to create a design matrix from a model formula and a data frame. It takes a formula object and a data frame as arguments and returns a design matrix based on the formula. The design matrix includes the predictor variables specified in the formula, as well as any transformations or interactions specified in the formula.

# The model.matrix() function is particularly useful for creating a design matrix for regression models that involve categorical variables. It automatically handles the encoding of categorical variables into appropriate dummy variables or contrasts based on the coding scheme specified in the formula. The resulting design matrix can then be used as input to various regression models in R.

# End ------------

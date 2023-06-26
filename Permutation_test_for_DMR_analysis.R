# Permutation test to get pvalues for DMRs:

# Start ------------

# Loading libraries ------------

library(bsseq)
library(data.table)
library(limma)
library(parallel)

# Uploading BSseq object for all KO vs all control design

BSseq.obj <- readRDS("/media/ko/New Volume/Documents/Nanopore_data/modbam2bed_outputs/UCSC_mouse_genome_based_analysis/mapping_quality_30_and_400bp_length_filtering_results/DMR_results_all_KOs_versus_all_control_groups/Last_analysis_with_last_pipeline/merged_BSseq_smoothed_filtered.rds")

# Creating design matrix ------------

pData(BSseq.obj)$metadata_2[,1] # This returns annotation of samples, that is, Control and KOs!
# metada_2 annotation is a data frame but we need to have character strings instead of data frame

BSseq.obj$group <- pData(BSseq.obj)$metadata_2[,1] # create "group" column in BSseq object desribing samples
# as control and knockouts

design <- model.matrix(~ group, colData(BSseq.obj))
colnames(design) <- gsub("group", "", colnames(design))
contrasts <- diag(rep(1, ncol(design)))[, -1, drop = FALSE]
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
  maxGap.sd = 10 ^ 8,
  maxGap.dmr = 300,
  type = "dmrs",
  # NOTE: This chunksize is a good choice if your data are chunked row-wise
  #       (e.g., 100 Mb blocks).
  #chunksize = chunkdim(assay(BSseq.obj, "coef"))[1],
  mc.cores = as.integer(Sys.getenv("NSLOTS")))

# Save results as rds file ------------

saveRDS(
  object = fstat_pipeline$dmrs[fstat_pipeline$dmrs$fwer < 50, ],
  file = here("DMRs", "objects", "general_CG-DMRs.Phase1.rds"))

# End ------------
# Permutation test to get pvalues for DMRs:

# Start ------------

# Loading libraries ------------

library(bsseq)
library(data.table)
library(limma)
library(parallel)
library(permute)

# Uploading BSseq object for all KO vs all control design

BSseq.obj <- readRDS("/hpcdata/Mimir/shared/Juan_ONT_data/processed_data/results/allKOs_vs_all_controls/BSseq_obj/merged_BSseq_smoothed_filtered.rds")

# Creating design matrix ------------

pData(BSseq.obj)$metadata_2[,1] # This returns annotation of samples, that is, Control and KOs!
# metada_2 annotation is a data frame but we need to have character strings instead of data frame

BSseq.obj$group <- pData(BSseq.obj)$metadata_2[,1] # create "group" column in BSseq object describing samples
# as control and knockouts

# Create design matrix to assing control samples "1", where knockout samples 0.
# Design matrix is kind of descriptive data for linear model.

design <- model.matrix(~ group, colData(BSseq.obj))
colnames(design) <- gsub("group", "", colnames(design))

# Create Constrasts matrix

contrasts <- diag(rep(1, ncol(design)))
rownames(contrasts) <- colnames(design)

# Create all functions needed in permutation test start ----- (all custom scripts obtained from "https://github.com/hansenlab/bsseq/blob/devel/R/permutations.R")

# Create custom function getNullDistribution_BSmooth.tstat()

getNullDistribution_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2,
                                              estimate.var, local.correct,
                                              cutoff, stat, maxGap, mc.cores = 1) {
    stopifnot(nrow(idxMatrix1) == nrow(idxMatrix2))
    message(sprintf("[getNullDistribution_BSmooth.tstat] performing %d permutations\n", nrow(idxMatrix1)))
    nullDist <- mclapply(1:nrow(idxMatrix1), function(ii) {
        ptime1 <- proc.time()
        BS.tstat <- BSmooth.tstat(BSseq, estimate.var = estimate.var,
                                  group1 = idxMatrix1[ii,],
                                  group2 = idxMatrix2[ii,],
                                  local.correct = local.correct, maxGap = 10^8,
                                  verbose = FALSE)
        dmrs0 <- dmrFinder(BS.tstat, stat = stat, cutoff = cutoff, maxGap = maxGap)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        message(sprintf("[getNullDistribution_BSmooth.tstat] completing permutation %d in %.1f sec\n", ii, stime))
        dmrs0
    }, mc.cores = min(nrow(idxMatrix1), mc.cores), mc.preschedule = FALSE)
    nullDist
}

# Create custom function getNullDistribution_BSmooth.fstat()

getNullDistribution_BSmooth.fstat <- function(BSseq,
                                              idxMatrix,
                                              design, contrasts,
                                              coef = NULL, cutoff,
                                              maxGap.sd, maxGap.dmr,
                                              mc.cores = 1) {

    message(sprintf("[getNullDistribution_BSmooth.fstat] performing %d permutations\n",
                    nrow(idxMatrix)))
    # NOTE: Using mc.preschedule = TRUE
    # TODO: Need some protection for when a core(s) error, otherwise nullDist
    #       contains a mixture of data frame, NULL, and try-error objects
    #       (which will subesequently break when passed to getFWER.fstat())
    nullDist <- mclapply(seq_len(nrow(idxMatrix)), function(ii) {
        ptime1 <- proc.time()
        # NOTE: More efficient to permute design matrix using idxMatrix[ii, ] than
        #       permute the raw data with BSseq[, idxMatrix[ii, ]]
        bstat <- BSmooth.fstat(BSseq, design = design[idxMatrix[ii, ], ],
                               contrasts = contrasts)
        bstat <- smoothSds(bstat, maxGap = maxGap.sd)
        bstat <- computeStat(bstat, coef = coef)
        dmrs0 <- dmrFinder(bstat, cutoff = cutoff, maxGap = maxGap.dmr)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        message(sprintf("[getNullDistribution_BSmooth.fstat] completing permutation %d in %.1f sec\n", ii, stime))
        dmrs0
    }, mc.cores = min(nrow(idxMatrix), mc.cores), mc.preschedule = TRUE)
    nullDist
}

# Create custom function permuteAll()

permuteAll <- function(nperm, design) {
    message(sprintf("[permuteAll] performing %d unrestricted permutations of the design matrix\n", nperm))

    CTRL <- how(nperm = nperm)
    # NOTE: shuffleSet() returns a nperm * nrow(design) matrix of permutations
    idxMatrix <- shuffleSet(n = nrow(design), control = CTRL)
}

# Create custom function getNullBlocks_BSmooth.tstat()

getNullBlocks_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                          mc.cores = 1) {
    getNullDistribution_BSmooth.tstat(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = FALSE,
                       estimate.var = estimate.var,
                       cutoff = c(-2,2), stat = "tstat", maxGap = 10000,
                       mc.cores = mc.cores)
}

# Create custom fuction getNullDmrs()

getNullDmrs_BSmooth.tstat <- function(BSseq, idxMatrix1, idxMatrix2, estimate.var = "same",
                        mc.cores = 1) {
    getNullDistribution_BSmooth.tstat(BSseq = BSseq, idxMatrix1 = idxMatrix1,
                       idxMatrix2 = idxMatrix2, local.correct = TRUE,
                       estimate.var = estimate.var,
                       cutoff = c(-4.6,4.6), stat = "tstat.corrected", maxGap = 300,
                       mc.cores = mc.cores)
}


# Create custom function subsetDmrs()

subsetDmrs <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    out <- xx[ xx[,"n"] >= 3 & abs(xx[, "meanDiff"]) > 0.1 &
                  xx[, "invdensity"] <= 300, ]
    if(nrow(out) == 0)
        return(NULL)
    out
}

# Create custom function subsetBlocks()

subsetBlocks <- function(xx) {
    if(is.null(xx) || is(xx, "try-error"))
        return(NULL)
    out <- subset(xx, width >= 10000)
    if(nrow(out) == 0)
        return(NULL)
    out
}

# Create custom function getFWER()

getFWER <- function(null, type = "blocks") {
    reference <- null[[1]]
    null <- null[-1]
    null <- null[!sapply(null, is.null)]
    better <- sapply(1:nrow(reference), function(ii) {
        # meanDiff <- abs(reference$meanDiff[ii])
        areaStat <- abs(reference$areaStat[ii])
        width <- reference$width[ii]
        n <- reference$n[ii]
        if (type == "blocks") {
            out <- sapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                        # nulldist$width >= width)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$width >= width)
            })
        }
        if (type == "dmrs") {
            out <- sapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                #     nulldist$n >= n)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$n >= n)
            })
        }
        sum(out)
    })
    better
}


# Create custom function getFWER.fstat()

# NOTE: Identical to getFWER() except uses areaStat rather than meanDiff
#       to compare regions.
getFWER.fstat <- function(null, type = "blocks") {
    reference <- null[[1]]
    null <- null[-1]
    null <- null[!sapply(null, is.null)]
    # TODO: Will break if null == list(), which can occur in practice (although
    #       rarely).
    better <- sapply(seq_len(nrow(reference)), function(ii) {
        # meanDiff <- abs(reference$meanDiff[ii])
        areaStat <- abs(reference$areaStat[ii])
        width <- reference$width[ii]
        n <- reference$n[ii]
        if (type == "blocks") {
            out <- vapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                # nulldist$width >= width)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$width >= width)
            }, logical(1L))
        }
        if (type == "dmrs") {
            out <- vapply(null, function(nulldist) {
                # any(abs(nulldist$meanDiff) >= meanDiff &
                #     nulldist$n >= n)
                any(abs(nulldist$areaStat) >= areaStat &
                        nulldist$n >= n)
            }, logical(1L))
        }
        sum(out)
    })
    better
}


# Create custom function makeIdxMatrix()

# TODO: Simplify makeIdxMatrix() by using permute package
makeIdxMatrix <- function(group1, group2, testIsSymmetric = TRUE, includeUnbalanced = TRUE) {
    groupBoth <- c(group1, group2)
    idxMatrix1 <- NULL
    subsetByMatrix <- function(vec, mat) {
        apply(mat, 2, function(xx) vec[xx])
    }
    combineMat <- function(mat1, mat2) {
        tmp <- lapply(1:nrow(mat1), function(ii) {
            t(apply(mat2, 1, function(xx) { c(mat1[ii,], xx) }))
        })
        do.call(rbind, tmp)
    }
    if(length(group1) == 1 && length(group1) == 1) {
        if(testIsSymmetric)
            idxMatrix1 <- as.matrix(group1)
        else
            idxMatrix1 <- as.matrix(c(group1, group2))
    }
    if(length(group1) == 2 && length(group2) == 2) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(matrix(group1[1], ncol = 1),
                                           matrix(group2, ncol = 1)))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(matrix(group1, ncol = 1), matrix(group2, ncol = 1)))
        }
    }
    if(length(group1) == 3 && length(group1) == 3) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           as.matrix(group2, ncol = 1)))

        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           as.matrix(group2, ncol = 1)),
                                combineMat(as.matrix(group1, ncol = 1),
                                           subsetByMatrix(group2, combinations(3,2))))
        }
    }
    if(length(group1) == 4 && length(group1) == 4) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(3,2)),
                                           subsetByMatrix(group2, combinations(4,2))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(4,2)),
                                           subsetByMatrix(group2, combinations(4,2))))
        }
        if(includeUnbalanced) {
            newMatrix <- combineMat(subsetByMatrix(group1, combinations(4,3)),
                                    as.matrix(group2, ncol = 1))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix)
        }
        if(includeUnbalanced && !testIsSymmetric) {
            newMatrix <- combineMat(as.matrix(group1, ncol = 1),
                                    subsetByMatrix(group2, combinations(4,3)))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix)
        }
    }
    if(length(group1) == 5 && length(group1) == 5) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(5, 3)),
                                           subsetByMatrix(group2, combinations(5, 2))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(5, 3)),
                                           subsetByMatrix(group2, combinations(5, 2))),
                                combineMat(subsetByMatrix(group1, combinations(5, 2)),
                                           subsetByMatrix(group2, combinations(5, 3))))
        }
        if(includeUnbalanced) {
            idxMatrix1 <- rbind(idxMatrix1,
                                combineMat(subsetByMatrix(group1, combinations(5,4)),
                                           as.matrix(group2, ncol = 1)))
        }
        if(includeUnbalanced && !testIsSymmetric) {
            idxMatrix1 <- rbind(idxMatrix1,
                                combineMat(as.matrix(group1, ncol = 1),
                                           subsetByMatrix(group2, combinations(5,4))))
        }
    }
    if(length(group1) == 6 && length(group1) == 6) {
        if(testIsSymmetric) {
            idxMatrix1 <- rbind(group1,
                                combineMat(subsetByMatrix(group1, combinations(5,3)),
                                           subsetByMatrix(group2, combinations(6,3))))
        } else {
            idxMatrix1 <- rbind(group1, group2,
                                combineMat(subsetByMatrix(group1, combinations(6,3)),
                                           subsetByMatrix(group2, combinations(6,3))))
        }
        if(includeUnbalanced) {
            newMatrix1 <- combineMat(subsetByMatrix(group1, combinations(6,4)),
                                    subsetByMatrix(group2, combinations(6,2)))
            newMatrix2 <- combineMat(subsetByMatrix(group1, combinations(6,5)),
                                    as.matrix(group2, ncol = 1))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix1, newMatrix2)
        }
        if(includeUnbalanced && !testIsSymmetric) {
            newMatrix1 <- combineMat(subsetByMatrix(group1, combinations(6,2)),
                                     subsetByMatrix(group2, combinations(6,4)))
            newMatrix2 <- combineMat(as.matrix(group1, ncol = 1),
                                     subsetByMatrix(group2, combinations(6,5)))
            idxMatrix1 <- rbind(idxMatrix1, newMatrix1, newMatrix2)
        }
    }
    if(is.null(idxMatrix1))
        stop("unable to handle this combination of 'group1', 'group2' and 'testIsSymmetric'")
    rownames(idxMatrix1) <- NULL
    idxMatrix2 <- do.call(rbind, lapply(1:nrow(idxMatrix1), function(ii) {
        setdiff(groupBoth, idxMatrix1[ii,])
    }))
    return(list(idxMatrix1 = idxMatrix1, idxMatrix2 = idxMatrix2))
}

# Create custom function fstat.pipeline()

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

# Create custom function BSmooth.fstat()

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

# Create custom function smoothSds() ------------

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

# Create custom makeClusters() function ------------

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

# Create custom computeStat() function ------------

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

# Create all functions needed in permutation test end -----

# Do permutation test start ------------

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

# Do permutation test end ------------

# Save results either

saveRDS(object = fstat_pipeline$dmrs, file = "/hpcdata/Mimir/shared/Juan_ONT_data/processed_data/results/allKOs_vs_all_controls/Permutation_test_to_get_pvalues_for_DMRs/allKO_vs_allcontrols_all_DMRs_permutation_test_results.rds")

# or

write.csv(subset(fstat_pipeline$dmrs, n >= 3 & abs(maxDiff) >= 0.1), "/hpcdata/Mimir/shared/Juan_ONT_data/processed_data/results/allKOs_vs_all_controls/Permutation_test_to_get_pvalues_for_DMRs/allKO_vs_allcontrols_permutation_test_DMR_results.csv", quote = F, row.names = F)

# End ------------

# Notes from Kasper:

# fstat has no local.correct and estimate.var arguments so its not same with tstat for DMR analysis. When you do DMR analysis with fstat, it would give different DMRs than tsats gives.
# Put k argument and local.correct argument to getNullDistribution_BSmooth.tstat() function for permutation test.
# Go permuteAll() function to see how we can reshuffle or resampling samples (knockout and control samples) for permutation.
# WGBS analysis, this permutation pipeline is okay with equal number of samples, but we have different number of samples (13 knockout and 5 controls) so we can't generate design or idx matrix for permutation.
# We need to create new function for permutation addressing asymmetric permutation where number of samples is different.
# We function (tstat and fstat) should be merged somehow for permutation.
# Each row in idxmatrix is permutation of samples.
# getFWER function for calculation of pvalues. That's why Kasper decided to put <50 cutoff for it because we're permutating our data 1000 times to get pvalues lower than 0.05 and when we put <50 FWER cutoff which means
# 50/1000 <- 0.05.
# We need to permute samples ids!
# Check shuffleset() and permuteAll() function for permutation of sample ids.
                  

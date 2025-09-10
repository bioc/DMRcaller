#' Compute Co-Methylation Between CpG Sites Within regions
#'
#' @title Compute pairwise co-methylation statistics for cytosine sites within regions
#' @description
#' \code{computeCoMethylation()} calculates pairwise co-methylation between all cytosine sites
#' within each given region, using ONT methylation calls annotated to each site.
#' For each pair of cytosines within the same strand and PMD, it builds a 2x2 contingency table
#' reflecting the overlap state of reads (both methylated, only one methylated, or neither),
#' performs a statistical test (Fisher's exact by default), and reports FDR-adjusted p-values.
#'
#' @param methylationData A \code{GRanges} object containing cytosine sites, annotated with
#'   per-site ONT methylation calls (columns \code{ONT_Cm}, \code{ONT_C}, \code{readsN}, etc).
#' @param regions A \code{GRanges} object with list including genomic context such as gene and/or transposable
#' elements coordinates which possibly have DMRs, VMRs or PMDs.
#' @param minDistance Minimum distance (in bp) between two cytosines to consider for co-methylation (default: 150).
#' @param maxDistance Maximum distance (in bp) between two cytosines to consider (default: 1000).
#' @param minCoverage Minimum read coverage required for both cytosines in a pair (default: 4).
#' @param pValueThreshold FDR-adjusted p-value threshold for reporting significant co-methylation (default: 0.01).
#' @param test Statistical test to use for co-methylation (\code{"fisher"}
#' for Fisher's exact [default], or \code{"permutation"} for chi-squared).
#' NOTE: highly recommended to do parallel when use permutation test.
#' @param alternative indicates the alternative hypothesis and must be one of
#' \code{"two.sided"}, \code{"greater"} or \code{"less"}. You can specify just
#' the initial letter. Only used in the 2 by 2 case. This is used only for
#' Fisher's test.
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#' This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @return A list of length equal to \code{regions}, where each entry is a \code{GInteractions} object
#'   of significant cytosine pairs (by strand), annotated with:
#'   \describe{
#'     \item{C1_C2}{number of reads methylated at both cytosines}
#'     \item{C1_only}{number methylated at only first cytosine}
#'     \item{C2_only}{number methylated at only second cytosine}
#'     \item{neither}{number methylated at neither cytosines}
#'     \item{strand}{The DNA strand ("+" or "-") on which the two CpGs reside.}
#'     \item{genomic_position}{The original region (from \code{regions}) containing
#'     this cytosines pair, formatted in UCSC or IGV style, e.g. \code{"chr1:1522971-1523970"}.}
#'     \item{p.value}{FDR-adjusted p-value for co-methylation association}
#'   }
#'
#' @details
#' Pairwise tests are performed separately for each strand (+ and -) within each region.
#' FDR correction is performed for all pairs within each region and strand.
#'
#' @seealso \code{\link{readONTbam}}, \code{\link{computePMDs}},
#' \code{\link{ontSampleGRangesList}}
#'
#' @examples
#' \dontrun{
#' # load the ONT methylation data and PMD data
#' data(ont_gr_GM18870_chr1_PMD_bins_1k)
#' data(ont_gr_GM18870_chr1_sorted_bins_1k)
#'
#' # compute the co-methylations with Fisher's exact test
#' coMetylationFisher <- computeCoMethylation(
#'   ont_gr_GM18870_chr1_sorted_bins_1k,
#'   regions = ont_gr_GM18870_chr1_PMD_bins_1k[1:4],
#'   minDistance = 150,
#'   maxDistance = 1000,
#'   minCoverage = 4,
#'   pValueThreshold = 0.01,
#'   test = "fisher",
#'   parallel = FALSE)
#'
#' # compute the co-methylations with Permuation test
#' coMetylationPermutation <- computeCoMethylation(
#'   ont_gr_GM18870_chr1_sorted_bins_1k,
#'   regions = ont_gr_GM18870_chr1_PMD_bins_1k[1:4],
#'   minDistance = 150,
#'   maxDistance = 1000,
#'   minCoverage = 4,
#'   pValueThreshold = 0.01,
#'   test = "permutation",
#'   parallel = FALSE) # highly recommended to set as TRUE
#' }
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import InteractionSet
#' @import BiocParallel
#' @export
computeCoMethylation <- function(methylationData,
                                          regions,
                                          minDistance = 150,
                                          maxDistance = 1000,
                                          minCoverage = 4,
                                          pValueThreshold = 0.01,
                                          alternative = "two.sided",
                                          test = "fisher",
                                          parallel = FALSE,
                                          BPPARAM  = NULL){
  ##Parameters checking
  cat("Parameters checking ...\n")

  # generate the BPPARAM value if set as parallel
  if (parallel == TRUE){
    BPPARAM <- suppressWarnings(.validateBPPARAM(BPPARAM, progressbar = TRUE))
    cat("Current parallel setting, BPPARAM: ",
        capture.output(BPPARAM),sep = "\n")
  }else{
    # Force serial execution
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
    cat("Current parallel setting, BPPARAM: ",
        capture.output(BPPARAM),sep = "\n")
  }

  .validateMethylationData(methylationData, variableName="methylationData")
  regions <- .validateGRanges(regions, methylationData)

  .validateCoMethylationStatTest(test)

  .validateCoMethylationAlternative(alternative)

  .stopIfNotAll(c(!is.null(minDistance), is.numeric(minDistance), minDistance >= 0),
                " minDistance needs to be a numeric value higher or equal to 0")

  .stopIfNotAll(c(!is.null(maxDistance), is.numeric(maxDistance), maxDistance >= 0),
                " maxDistance needs to be a numeric value higher or equal to 0")

  .stopIfNotAll(c(maxDistance >= minDistance),
                " maxMethylation should be higher than minMethylation value")

  .stopIfNotAll(c(.isInteger(minCoverage, positive=TRUE)),
                " the minimum gap between PMDs is an integer higher or equal to 0")

  .stopIfNotAll(c(!is.null(pValueThreshold),
                  is.numeric(pValueThreshold),
                  pValueThreshold > 0,
                  pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")

  results_gil =  vector(mode = "list", length = length(regions))
  regions <- .generateGRangesName(regions)
  t_start <- Sys.time()
  n_total <- length(regions)
  for (i in seq_along(regions)) {
    elapsed <- difftime(Sys.time(), t_start, units = "secs")
    cat(sprintf("[computeCoMethylation] PMD %d/%d (%.1f%%) at %s | Elapsed: %.1f sec\n",
                i, n_total, 100*i/n_total, format(Sys.time(), "%H:%M:%S"), as.numeric(elapsed)))
    bin <- regions[i]
    in_bin <- methylationData[queryHits(findOverlaps(methylationData, bin))]
    if (length(in_bin) < 2) next
    all_strands <- unique(as.character(strand(in_bin)))
    bin_results <- list()
    for (s in all_strands) {
      in_bin_strand <- in_bin[strand(in_bin) == s]
      bin_results[[s]] <- .processStrandPairs(
        in_bin_strand,
        minDistance = minDistance,
        maxDistance = maxDistance,
        minCoverage = minCoverage,
        pValueThreshold = pValueThreshold,
        test = test,
        alternative = alternative
      )
    }
    merged <- do.call(c, bin_results)
    results_gil[[i]] <- merged
  }
  names(results_gil) <- regions$genomic_position

  # keep only those bins where we got a valid GInteractions back
  keep <- vapply(results_gil,
                 function(x) !is.null(x) & length(x) > 0L,
                 logical(1))
  results_gil <- results_gil[keep]

  results_gi <- do.call(c, Map(function(strand_list, region_name) {
    # 1) pull out the two strand‐specific GInteractions
    gi_pos <- strand_list[["+"]]
    gi_neg <- strand_list[["-"]]
    # guard against NULL
    if (is.null(gi_pos)) gi_pos <- GInteractions()
    if (is.null(gi_neg)) gi_neg <- GInteractions()

    mcols(gi_pos)$strand <- rep("+",length(gi_pos))
    mcols(gi_neg)$strand <- rep("-",length(gi_neg))
    # 2) concatenate them
    gi <- c(gi_pos, gi_neg)
    # 3) add a column ‘genomic_position’ with the region’s name
    mcols(gi)$genomic_position <- rep(region_name, length(gi))
    gi
  }, results_gil, names(results_gil)))

  result <- Reduce(c, results_gi)

  # Adjust the p.value using the FDR method
  mcols(result)$p.value <- stats::p.adjust(result$pre_p.value, method = 'fdr')
  mcols(result)$pre_p.value <- NULL
  keep_final <- result$p.value < pValueThreshold
  result <- result[which(result$p.value < pValueThreshold & !is.na(result$p.value))]


  total_elapsed <- difftime(Sys.time(), t_start, units = "secs")
  cat(sprintf("[computeCoMethylation] Done! Total elapsed time: %.1f sec\n", as.numeric(total_elapsed)))
  return(result)
}

# --------------
#  Helper: compute the co-methylation between the cpg in regions
# --------------
# 1. Calculate contingency table for two CpGs
.getCoMethylationContingencyTable <- function(c1, c2) {
  split_ids <- function(x) {
    v <- as.character(unlist(x))
    if (length(v) == 0L) return(character(0))
    v <- v[!is.na(v) & nzchar(v)]
    if (length(v) == 0L) return(character(0))
    out <- unlist(strsplit(v, ",", fixed = TRUE), use.names = FALSE)
    unique(out[nzchar(out)])
  }

  # Per-CpG sets: methylated (M) and unmethylated (U) read IDs
  M1 <- split_ids(c1$ONT_Cm)
  U1 <- split_ids(c1$ONT_C)
  M2 <- split_ids(c2$ONT_Cm)
  U2 <- split_ids(c2$ONT_C)

  R <- intersect(union(M1, U1), union(M2, U2))

  C1_C2   <- length(intersect(R, intersect(M1, M2)))
  C1_only <- length(intersect(R, intersect(M1, U2)))
  C2_only <- length(intersect(R, intersect(U1, M2)))
  neither <- length(intersect(R, intersect(U1, U2)))

  list(C1_C2 = C1_C2, C1_only = C1_only, C2_only = C2_only, neither = neither)

}

.getMethylationStatusVectors <- function(c1, c2) {
  reads1_meth <- unlist(strsplit(unlist(c1$ONT_Cm), split = ","))
  reads1_unmeth <- unlist(strsplit(unlist(c1$ONT_C), split = ","))
  reads2_meth <- unlist(strsplit(unlist(c2$ONT_Cm), split = ","))
  reads2_unmeth <- unlist(strsplit(unlist(c2$ONT_C), split = ","))
  all_reads <- unique(c(reads1_meth, reads1_unmeth, reads2_meth, reads2_unmeth))
  status1 <- ifelse(all_reads %in% reads1_meth, 1,
                    ifelse(all_reads %in% reads1_unmeth, 0, NA))
  status2 <- ifelse(all_reads %in% reads2_meth, 1,
                    ifelse(all_reads %in% reads2_unmeth, 0, NA))
  keep <- !is.na(status1) & !is.na(status2)
  status1_filt <- status1[keep]
  status2_filt <- status2[keep]
  return(list(status1 = status1_filt, status2 = status2_filt))
}

# 2. Statistical test block
.coMethylationFisherTest <- function(c1_c2, c1_only, c2_only,
                                     neither, alternative= "two.sided"){
    table <- matrix(
      c(c1_c2, c1_only, c2_only, neither),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("c1_yes", "c1_no"), c("c2_yes", "c2_no"))
    )
    return(stats::fisher.test(table, alternative = alternative)$p.value)
}

### score test need to change to other methods as detecting the discordant methylation
.coMethylationScoreTest <- function(c1_c2, c1_only, c2_only, neither){
    n_total <- c1_c2 + c1_only + c2_only + neither
    m1 <- c1_c2
    n1 <- c1_c2 + c1_only
    m2 <- c1_c2
    n2 <- c1_c2 + c2_only
    return(.scoreTest(m1, n1, m2, n2))
}

.coMethylationBinomTest <- function(c1_c2, c1_only, c2_only, neither, alternative= "two.sided"){
  n_total <- c1_c2 + c1_only + c2_only + neither
  concordance <- c1_c2 + neither
  return(stats::binom.test(concordance,n_total,p = 0.5, alternative = alternative)$p.value)
}

.coMethylationPermutationTest <- function(
    status1, status2,
    nperm = 1000,
    statfun = function(tab) stats::chisq.test(tab)$statistic){
  # Create the observed contingency table from the two status vectors
  obs_tab <- table(status1, status2)
  if (!all(dim(obs_tab) == c(2,2))) {
    return(NA_real_)
  }
  obs_stat <- suppressWarnings(statfun(obs_tab))

  # Permute status2 nperm times and recalculate the test statistic for each permutation
  perm_stats <- replicate(nperm, {
    perm_status2 <- sample(status2)
    perm_tab <- table(status1, perm_status2)
    if (!all(dim(perm_tab) == c(2,2))) return(0)
    suppressWarnings(statfun(perm_tab))
  })
  perm_stats <- perm_stats[!is.na(perm_stats)]
  pval <- (sum(perm_stats >= obs_stat) + 1) / (length(perm_stats) + 1)
  return(pval)
}

# 3. Process pairs in one strand; only output adjusted p-values in the results
.processStrandPairs <- function(in_bin_strand,
                                minDistance,
                                maxDistance,
                                minCoverage,
                                pValueThreshold,
                                test = "fisher",
                                alternative = "two.sided",
                                nperm = 1000,
                                BPPARAM = BiocParallel::bpparam()) {
  if (length(in_bin_strand) < 2) return(GInteractions())
  cpg_pos <- start(in_bin_strand)
  coverage <- in_bin_strand$readsN
  pair_idx <- combn(seq_along(in_bin_strand), 2)
  dist_bp <- abs(cpg_pos[pair_idx[1,]] - cpg_pos[pair_idx[2,]])
  keep <- which(
    dist_bp >= minDistance & dist_bp <= maxDistance &
      coverage[pair_idx[1,]] >= minCoverage & coverage[pair_idx[2,]] >= minCoverage
  )
  if (length(keep) == 0) return(GInteractions())

  idx1 <- pair_idx[1, keep]
  idx2 <- pair_idx[2, keep]
  n_pairs <- length(idx1)
  results <- GInteractions(
    granges(in_bin_strand[idx1]),
    granges(in_bin_strand[idx2])
  )

  # Updated meta columns
  C1_C2 <- integer(n_pairs)
  C1_only <- integer(n_pairs)
  C2_only <- integer(n_pairs)
  neither <- integer(n_pairs)
  pvals <- numeric(n_pairs)


  # Parallel computation of contingency table statistics for all CpG pairs
  pair_stats <- bplapply(seq_len(n_pairs), function(i) {
    cont <- .getCoMethylationContingencyTable(in_bin_strand[idx1[i]], in_bin_strand[idx2[i]])
    if (test == "permutation") {
      metStatus <- .getMethylationStatusVectors(in_bin_strand[idx1[i]], in_bin_strand[idx2[i]])
      pval <- .coMethylationPermutationTest(metStatus$status1, metStatus$status2, nperm = nperm)
    } else if (test == "fisher"){
      pval <- .coMethylationFisherTest(cont$C1_C2, cont$C1_only, cont$C2_only, cont$neither, alternative)
    } else if (test == "score"){
      pval <- .coMethylationScoreTest(cont$C1_C2, cont$C1_only, cont$C2_only, cont$neither)
    } else if (test == "binom"){
      pval <- .coMethylationBinomTest(cont$C1_C2, cont$C1_only, cont$C2_only, cont$neither, alternative)
    }
    c(C1_C2 = cont$C1_C2, C1_only = cont$C1_only, C2_only = cont$C2_only, neither = cont$neither, pval = pval)
  }, BPPARAM = BPPARAM)

  pair_stats_mat <- do.call(rbind, pair_stats)
  if (is.null(dim(pair_stats_mat))) pair_stats_mat <- t(as.matrix(pair_stats_mat))

  # Extract each statistic as a numeric vector for all pairs
  C1_C2    <- pair_stats_mat[, "C1_C2"]
  C1_only  <- pair_stats_mat[, "C1_only"]
  C2_only  <- pair_stats_mat[, "C2_only"]
  neither  <- pair_stats_mat[, "neither"]
  pvals    <- pair_stats_mat[, "pval"]

  # Subset the results and assign metadata columns to GInteractions object
  results$pre_p.value  <- pvals
  results$C1_C2    <- C1_C2
  results$C1_only  <- C1_only
  results$C2_only  <- C2_only
  results$neither  <- neither

  return(results)
}

### generate the genomic_position formatted by the UCSC or IGV style
.generateGRangesName <- function(regions) {
  mcols(regions)$genomic_position <- paste0(seqnames(regions), ":",start(regions), "-", end(regions))
  return(regions)
}

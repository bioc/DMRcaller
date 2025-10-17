#' Compute Co-Methylated Regions (CMRs)
#'
#' @title Compute pairwise co-methylation statistics between regions
#' @description
#' \code{computeCoMethylatedRegions()} calculates pairwise correlation statistics for methylation levels
#' across defined genomic regions (e.g., PMDs, Enhancer binding sites). For each region pair within the
#' specified distance range, the function computes per-read methylation proportions and
#' performs correlation testing (Pearson, Spearman, or Kendall). Pairs with strong correlations
#' (beyond user-defined thresholds) and significant p-values (FDR-adjusted) are returned.
#'
#' @param methylationData A \code{GRanges} object containing cytosine sites, annotated with
#' ONT methylation metadata (columns \code{ONT_Cm}, \code{ONT_C}, etc.).
#' @param regions A \code{GRanges} object defining genomic regions (e.g., PMDs, Enhancer binding sites) to
#' evaluate for CMRs.
#' @param minDistance Minimum genomic distance (in base pairs) between two regions to be considered (default: 500).
#' @param maxDistance Maximum genomic distance (in base pairs) between two regions (default: 50,000).
#' @param minCoverage Minimum number of shared reads (per region pair) required to compute correlation (default: 4).
#' @param pValueThreshold Significance threshold for FDR-adjusted p-values (default: 0.05).
#' @param correlation_test Statistical method to compute correlation; must be one of \code{"pearson"}, \code{"spearman"}, or \code{"kendall"} (default: \code{"pearson"}).
#' @param minCorrelation Minimum allowed correlation value for a significant result (must be in between -1 and 0; default: -0.5).
#' @param maxCorrelation Maximum allowed correlation value for a significant result (must be in between 0 and 1; default: 0.5).
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#' This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' 
#' @return A \code{GInteractions} object, where each row represents a significantly correlated pair of genomic regions
#' from the input \code{regions}. The anchors of each interaction correspond to original regions,
#' and their genomic coordinates are retained in the \code{anchor1} and \code{anchor2} slots.
#'
#' Additionally, a \code{genomic_position} meta-column is included to indicate the original coordinate ranges
#' (in UCSC/IGV format) for each interaction, aiding downstream interpretation or visualisation.
#'
#' Each interaction is annotated with:
#' \describe{
#'   \item{correlation}{Correlation coefficient between the two regions}
#'   \item{coverage}{Number of shared reads used for correlation}
#'   \item{p.value}{FDR-adjusted p-value for the correlation test}
#' }
#' @details
#' The function first identifies all region pairs within the user-defined distance range. For each pair,
#' it calculates methylation proportions per read across both regions, extracts common reads,
#' and tests correlation using the selected method. FDR correction is applied globally across all region pairs.
#'
#' @seealso \code{\link{readONTbam}}, \code{\link{computePMDs}},
#' \code{\link{ontSampleGRangesList}}
#'
#' @examples
#' # Load methylation data and PMD regions
#' data("ont_gr_GM18870_chr1_sorted_bins_1k")
#' data("ont_gr_GM18870_chr1_PMD_bins_1k")
#' 
#' # Compute highly correlated regions (Pearson)
#' coMethylationRegion_pearson <- computeCoMethylatedRegions(ont_gr_GM18870_chr1_sorted_bins_1k,
#'                                                           ont_gr_GM18870_chr1_PMD_bins_1k[1:5],
#'                                                           minDistance = 500,
#'                                                           maxDistance = 50000,
#'                                                           minCoverage = 4,
#'                                                           pValueThreshold = 0.05,
#'                                                           correlation_test = "pearson",
#'                                                           minCorrelation = -0.5,
#'                                                           maxCorrelation = 0.5,
#'                                                           parallel = FALSE,
#'                                                           BPPARAM  = NULL)
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import InteractionSet
#' @import BiocParallel
#' @importFrom stats cor.test
#' @export

# correlation between the methylated regions
## Correlation between Methylated Regions
computeCoMethylatedRegions <- function(methylationData,
                                       regions,
                                       minDistance = 500,
                                       maxDistance = 50000,
                                       minCoverage = 4,
                                       pValueThreshold = 0.05,
                                       correlation_test = "pearson",
                                       minCorrelation = -0.5,
                                       maxCorrelation = 0.5,
                                       parallel = FALSE,
                                       BPPARAM  = NULL,
                                       cores = NULL){
  ##Parameters checking
  cat("[computeCoMethylatedRegions] Parameters checking ...\n")
  
  # generate the BPPARAM value if set as parallel 
  if (parallel == TRUE){
    BPPARAM <- suppressWarnings(.validateBPPARAM(BPPARAM, progressbar = TRUE)) 
  }else{
    # Force serial execution
    BPPARAM <- BiocParallel::SerialParam(progressbar = TRUE)
  }
  # If cores argument is specified
  if (!is.null(cores)) {
    .stopIfNotAll(.isInteger(cores, positive = TRUE), 
                  "the number of cores used when computing the DMRs needs to be an integer higher or equal to 1.")
    
    # Check if user requested more cores than available
    if (cores > BPPARAM$workers) {
      warning(paste0("The number of requested cores (", cores, 
                     ") exceeds the available system cores (", BPPARAM$workers, 
                     "). Automatically setting cores to the maximum available (", 
                     BPPARAM$workers, ")."))
      cores <- BPPARAM$workers
    } else {
      message(paste0("Using user-specified core count: ", cores))
    }
    
    # Apply the final core number
    BPPARAM$workers <- cores
  } else {
    cores <- BPPARAM$workers
  }
  cat("Current parallel setting, BPPARAM: ", capture.output(BPPARAM),sep = "\n")
  
  .validateMethylationData(methylationData, variableName="methylationData")
  regions <- .validateGRanges(regions, methylationData)
  
  .validateCorrelationStatTest(correlation_test)
  
  .stopIfNotAll(c(!is.null(minDistance), is.numeric(minDistance), minDistance >= 0),
                " minDistance needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(maxDistance), is.numeric(maxDistance), maxDistance >= 0),
                " maxDistance needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(maxDistance >= minDistance),
                " maxDistance should be higher than minDistance value")
  
  .stopIfNotAll(c(!is.null(minCorrelation), is.numeric(minCorrelation), minCorrelation <= 0, minCorrelation >= -1),
                " minCorrelation needs to be a numeric value between -1 and 0")
  
  .stopIfNotAll(c(!is.null(maxCorrelation), is.numeric(maxCorrelation), maxCorrelation <= 1, maxCorrelation >= 0),
                " maxCorrelation needs to be a numeric value between 0 and 1")
  
  .stopIfNotAll(c(maxCorrelation >= minCorrelation),
                " maxCorrelation should be higher than minCorrelation value")
  
  .stopIfNotAll(c(.isInteger(minCoverage, positive=TRUE)),
                " the minimum gap between PMDs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(pValueThreshold),
                  is.numeric(pValueThreshold),
                  pValueThreshold > 0,
                  pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")
  
  t_start <- Sys.time()
  
  # Find the cases for comparison within minDistance and maxDistance
  cat("[computeCoMethylatedRegions] Count all case of pairs...\n")
  allPairs <- .makeRegionPairs(regions, minDistance=500, maxDistance=50000)

  n_total <- length(allPairs$distance)
  results_gi =  GInteractions()
  
  # Compute the correlation of proportion per read by pairs
  cat("[computeCoMethylatedRegions] Compute Correlated Methylated Regions...\n")
  results_list <- bplapply(seq_len(n_total), function(i) {
    bin1 <- regions[allPairs[i,1]]
    bin2 <- regions[allPairs[i,2]]
    
    comet_result <- .computeMoleculeCoMeth(
      methylationData,
      bin1,
      bin2,
      method = correlation_test, 
      minCoverage = minCoverage
    )
    
    if (!is.null(comet_result)) {
      result <- GInteractions(bin1, bin2)
      result$pre_p.value <- comet_result$pvalue
      result$coverage <- comet_result$n
      result$correlation <- comet_result$correlation
      return(result)
    } else {
      return(NULL)
    }
  }, BPPARAM = BPPARAM)
  
  results_gi <- do.call(c, results_list[!vapply(results_list, is.null, logical(1))])
  
  cat("[computeCoMethylatedRegions] Filter Correlated Methylated Regions \n")
  # Adjust the p.value using the FDR method
  results_gi$p.value <- stats::p.adjust(results_gi$pre_p.value, method = 'fdr')
  mcols(results_gi)$pre_p.value <- NULL
  
  results_gi <- results_gi[
    which(
      !is.na(results_gi$p.value) & results_gi$p.value < pValueThreshold &
      !is.na(results_gi$correlation) &
        (results_gi$correlation < minCorrelation |
        results_gi$correlation > maxCorrelation)
    )
  ]
  
  total_elapsed <- difftime(Sys.time(), t_start, units = "secs")
  cat(sprintf("[computeCoMethylatedRegions] Done! Total elapsed time: %.1f sec\n", as.numeric(total_elapsed)))
  return(results_gi)
}

# --------------
#  Helper: compute the correlation between the regions
# --------------
# 1. Find out which regions make pairs within specific distance
.makeRegionPairs <- function(regions, minDistance=500, maxDistance=50000) {
  # Generate all possible unique pairs of region indices (without repetition, no self-pairs)
  combs <- combn(seq_along(regions), 2, simplify=TRUE)  # matrix of 2 x N pairs
  
  # Compute distances for each pair
  dists <- apply(combs, 2, function(x) distance(regions[x[1]], regions[x[2]]))
  
  # Keep only valid pairs
  keep <- !is.na(dists) & dists >= minDistance & dists <= maxDistance
  
  # Return a data.frame with indices and distance
  data.frame(
    regionA = combs[1, keep],
    regionB = combs[2, keep],
    distance = dists[keep]
  )
}

# 2. Calculate Correlation from the two regions
.computeMoleculeCoMeth <- function(methylationData,
                                   regionA,
                                   regionB,
                                   method = "pearson", 
                                   minCoverage = 5){
  # subset reads overlapping region A and region B
  readsA <- methylationData[queryHits(findOverlaps(methylationData, regionA))]
  readsB <- methylationData[queryHits(findOverlaps(methylationData, regionB))]
  
  indexPropA <- .getProportionPerRead(readsA)
  indexPropB <- .getProportionPerRead(readsB)
  
  commonIds <- intersect(indexPropA$index_id,indexPropB$index_id)
  
  propCommonA =indexPropA$proportion[match(commonIds, indexPropA$index_id)]
  propCommonB =indexPropB$proportion[match(commonIds, indexPropB$index_id)]
  
  # compute correlation
  if(length(commonIds) >= minCoverage){
    cor_result <- cor.test(propCommonA,propCommonB, method = method)
    
      return(list(correlation=cor_result$estimate,
              pvalue=cor_result$p.value,
              n=length(commonIds)))
  }
}

# 3. Get the proportions and index from per-read
# Returns a list with:
# - index_id: character vector of unique read identifiers
# - proportion: numeric vector of methylation proportion per read (M / (M + U))
.getProportionPerRead <- function(methylationData){
  if (!("ONT_Cm" %in% names(mcols(methylationData))) |
      !("ONT_C"  %in% names(mcols(methylationData)))) {
    return(list(mean = NA_real_, sd = NA_real_, mean_w = NA_real_, sd_w = NA_real_))
  }
  
  # collect the sequence index from GRanges (ONT_Cm, ONT_C)
  split_ids <- function(x) {
    if (length(x) == 0L) return(character(0))
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(character(0))
    unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  }
  
  read_Cm_ids <- split_ids(unlist(methylationData$ONT_Cm))
  read_C_ids  <- split_ids(unlist(methylationData$ONT_C))
  total_idx   <- sort(unique(c(read_Cm_ids, read_C_ids)))
  total_idx   <- total_idx[nzchar(total_idx)]
  
  if (length(total_idx) == 0L) {
    return(list(mean = NA_real_, sd = NA_real_, mean_w = NA_real_, sd_w = NA_real_))
  }
  
  Cm_list <- lapply(mcols(methylationData)$ONT_Cm, split_ids)
  C_list  <- lapply(mcols(methylationData)$ONT_C,  split_ids)
  
  proportions_S <- numeric(length(total_idx))
  N_vec         <- numeric(length(total_idx))

  for (k in seq_along(total_idx)) {
    id <- total_idx[k]
    per_read_Cm_S <- vapply(Cm_list, function(vec) sum(vec == id), integer(1))
    per_read_C_S  <- vapply(C_list,  function(vec) sum(vec == id), integer(1))
    M <- sum(per_read_Cm_S)
    N <- M + sum(per_read_C_S)
    proportions_S[k] <- if (N > 0L) M / N else NA_real_
    N_vec[k]         <- N
  }

  # remove NA/NaN
  ok <- is.finite(proportions_S)
  total_idx <- total_idx[ok]
  p  <- proportions_S[ok]
  
  return(list(index_id = total_idx, proportion = p))
}

#' Checks whether the passed parameter is the statistical test for detecting Correlation in Regions
#'
#' @title Validate statistial test for Correlation between two Regions
#' @param test the statistical test used to call Correlation Between Methylated Regions (CMRs) (\code{"pearson"},
#'  \code{"kendall"} or \code{"spearman"} for Correlation test).
#' @keywords internal
#' @author Radu Zabet and Young Jun Kim
.validateCorrelationStatTest <- function(test){
  .stopIfNotAll(c(!is.null(test), is.character(test), length(test) == 1, test %in% c("pearson","kendall","spearman")),
                " Test needs to be one of the following \"pearson\", \"kendall\" or \"spearman\" for Correlation test" )
}

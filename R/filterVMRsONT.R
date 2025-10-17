#' Filter VMRs with ONT-specific variance tests and CI filters
#'
#' This function verifies whether a set of potential VMRs (e.g., genes,
#' transposons, CpG islands) are differentially methylated or not in ONT data,
#' adding per-read Wilcoxon and F-tests on per-site proportions, confidence interval
#' filtering, and optional variance-fold change cutoffs.
#'
#' @title Filter VMRs for ONT Data
#'
#' @param methylationData1 A \code{GRanges} of methylation calls for condition 1
#'   (see \code{\link{ontSampleGRangesList}}).
#' @param methylationData2 A \code{GRanges} of methylation calls for condition 2.
#' @param potentialVMRs A \code{GRanges} of candidate VMR regions (genes, TEs,
#'   CpG islands, etc.).
#' @param context Character string specifying cytosine context ("CG", "CHG",
#'   or "CHH").
#' @param pValueThreshold Numeric p-value threshold (0<value<1) for both
#'   Wilcoxon and F-tests after FDR adjustment.
#' @param minCytosinesCount Integer minimum number of cytosines per region.
#' @param minProportionDifference Numeric minimum methylation difference
#'   between conditions (0<value<1).
#' @param minReadsPerCytosine Integer minimum average coverage per cytosine.
#' @param ciExcludesOne Logical; if \code{TRUE}, filter out regions whose
#'   F-test 95\% confidence interval spans 1 (i.e., no significant variance change).
#' @param varRatioFc Optional; numeric fold-change cutoff on variance ratio
#'   (e.g., 2 for twofold variance difference). Regions with variance ratio
#'   outside \code{[1/varRatioFc, varRatioFc]} are kept when set.
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed \code{BPPARAM$workers}).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#'
#' @return A \code{GRanges} with the same ranges as \code{regions}, plus these metadata:
#' \describe{
#'   \item{sumReadsM1}{total methylated reads in condition 1}
#'   \item{sumReadsN1}{total reads in condition 1}
#'   \item{proportion1}{methylation proportion (sumReadsM1/sumReadsN1)}
#'   \item{variance1}{variance of per-read methylation proportions in condition 1}
#'   \item{sumReadsM2}{total methylated reads in condition 2}
#'   \item{sumReadsN2}{total reads in condition 2}
#'   \item{proportion2}{methylation proportion (sumReadsM2/sumReadsN2)}
#'   \item{variance2}{variance of per-read methylation proportions in condition 2}
#'   \item{cytosinesCount}{number of cytosines observed in each region}
#'   \item{wilcox_pvalue}{FDR adjusted p-value from Wilcoxon rank-sum test comparing per-read proportions}
#'   \item{f_pvalue}{FDR adjusted p-value from F-test comparing variances of per-read proportions}
#'   \item{var_ratio}{Ratio of variances (variance1 / variance2)}
#'   \item{wilcox_result}{Full \code{htest} object returned by \code{wilcox.test}}
#'   \item{F_test_result}{Full \code{htest} object returned by \code{var.test}}
#'   \item{direction}{a number indicating whether the region lost (-1)  or gain
#'     (+1) methylation in condition 2 compared to condition 1.}
#'   \item{regionType}{a string indicating whether the region lost (\code{"loss"})
#'     or gained (\code{"gain"}) methylation in condition 2 compared to condition 1.}
#'   \item{is_DMR}{logical; \code{TRUE} if region passed the \code{wilcox.test}}
#'   \item{is_VMR}{logical; \code{TRUE} if region passed the \code{var.test}}
#' }
#'
#' @details
#' For each potential VMR, per-site methylation proportions are aggregated per read,
#' then a two-sample Wilcoxon rank-sum test compares means (\code{wilcox_pvalue}), and
#' an F-test compares variances (\code{f_pvalue}). You may further filter by requiring
#' the 95% confidence interval from the F-test to exclude 1 (\code{ciExcludesOne}) and/or
#' apply a fold-change cutoff on the variance ratio (\code{varRatioFc}).
#'
#' @seealso \code{\link{readONTbam}},
#' \code{\link{computePMDs}}, \code{\link{computeCoMethylatedPositions}},
#' \code{\link{ontSampleGRangesList}}, \code{\link{GEs_hg38}}
#'
#' @examples 
#' \dontrun{
#' # load the ONT methylation data
#' data(ontSampleGRangesList)
#' # load the gene annotation data
#' data(GEs_hg38)
#'
#' # select the transcript
#' transcript <- GEs_hg38[which(GEs_hg38$type == "transcript")]
#'
#' # the regions where to compute the PMDs
#' regions <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' transcript <- transcript[overlapsAny(transcript, regions)]
#'
#' # filter genes that are differntially methylated in the two conditions
#' VMRsGenesCG <- filterVMRsONT(ontSampleGRangesList[["GM18501"]],
#'                ontSampleGRangesList[["GM18876"]], potentialVMRs = transcript,
#'                context = "CG", pValueThreshold = 0.01,
#'                minCytosinesCount = 4, minProportionDifference = 0.01,
#'                minReadsPerCytosine = 3, ciExcludesOne = TRUE,
#'                varRatioFc = NULL, parallel = TRUE) # parallel recommended
#' }
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @import stringr
#' @export

filterVMRsONT <- function(methylationData1,
                       methylationData2,
                       potentialVMRs,
                       context = "CG",
                       pValueThreshold = 0.01,
                       minCytosinesCount = 4,
                       minProportionDifference = 0.4,
                       minReadsPerCytosine = 3,
                       ciExcludesOne = TRUE,
                       varRatioFc = NULL, ## set as 2 if you want 2 fold-change varience cut-off
                       parallel = FALSE,
                       BPPARAM = NULL,
                       cores = NULL) {
  ##Parameters checking
  cat("Parameters checking ...\n")

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


  .validateMethylationData(methylationData1, variableName="methylationData1")
  .validateMethylationData(methylationData2, variableName="methylationData2")

  regions <- union(getWholeChromosomes(methylationData1),
                   getWholeChromosomes(methylationData2))

  .validateContext(context)

  .validateGRanges(potentialVMRs, generateGenomeWide=FALSE, variableName="potentialVMRs", minLength=NULL)

  .stopIfNotAll(c(!is.null(pValueThreshold), is.numeric(pValueThreshold), pValueThreshold > 0, pValueThreshold < 1),
                " the p-value threshold needs to be in the interval (0,1)")

  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)),
                " the minCytosinesCount is an integer higher or equal to 0")

  .stopIfNotAll(c(!is.null(minProportionDifference), is.numeric(minProportionDifference), minProportionDifference > 0, minProportionDifference < 1),
                " the minimum difference in methylation needs to be in the interval (0,1)")

  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)),
                " the minimum number of reads in a bin is an integer higher or equal to 0")

  .stopIfNotAll(is.logical(ciExcludesOne),
                " choose ciExcludesOne as TRUE if want to filter out regions whose F-test 95% confidence interval spans 1")

  if (!is.null(varRatioFc)) {
    .stopIfNotAll(c(is.numeric(varRatioFc), length(varRatioFc) == 1, varRatioFc >= 1),
    "`varRatioFc` must be a single numeric value >= 1 (e.g. 2 for a 2 times variance cutoff)")
  }


  regions <- reduce(regions)

  if(length(potentialVMRs) > 0){

    # extract the methylation data in the correct context
    cat("Extract methylation in the corresponding context \n")

    contextMethylationData1 <- methylationData1[methylationData1$context%in%context]
    rm(methylationData1)
    localContextMethylationData1 <- contextMethylationData1[queryHits(findOverlaps(contextMethylationData1, regions))]
    rm(contextMethylationData1)

    contextMethylationData2 <- methylationData2[methylationData2$context%in%context]
    rm(methylationData2)
    localContextMethylationData2 <- contextMethylationData2[queryHits(findOverlaps(contextMethylationData2, regions))]
    rm(contextMethylationData2)

    localContextMethylationData <- .joinMethylationData(localContextMethylationData1, localContextMethylationData2)
    rm(localContextMethylationData1, localContextMethylationData2)


    regionsList <- .splitGRangesEqualy(regions, cores)

    .filterVMRsLoop = function(i){
      computedVMRs <- GRanges()
      for(index in 1:length(regionsList[[i]])){
        currentRegion <- regionsList[[i]][index]


        cat("Computing VMRs at ",.printGenomicRanges(currentRegion),"\n")

        cat("Selecting data...\n")

        # Select the points in methylationData that we're interested in. These are the
        # points that lie within 'regions', as well as any that lie within
        # window.size of them.

        overlapsPotentialVMRs <- findOverlaps(potentialVMRs, currentRegion)
        if(length(overlapsPotentialVMRs) > 0){
          potentialVMRsLocal <- potentialVMRs[queryHits(overlapsPotentialVMRs)]

          localMethylationData <- localContextMethylationData[queryHits(findOverlaps(localContextMethylationData, currentRegion))]
          potentialVMRsLocal <- .analyseReadsInsideRegionsVMR(localMethylationData, potentialVMRsLocal)

          if(length(computedVMRs) == 0){
            computedVMRs <- potentialVMRsLocal
          } else{
            computedVMRs <- c(computedVMRs,potentialVMRsLocal)
          }
        }
      }
      return(computedVMRs)
    }

    # compute the VMRs
    if(cores > 1){
      cat("Compute the VMRs using ", cores, "cores\n")
      computedVMRs <- BiocParallel::bplapply(1:length(regionsList), .filterVMRsLoop, BPPARAM = BPPARAM)
    } else {
      computedVMRs <- lapply(1:length(regionsList), .filterVMRsLoop)
    }

    computedVMRs <-  unlist(GRangesList(computedVMRs))

    if(length(computedVMRs) > 0){
      computedVMRs$regionType <- rep("loss", length(computedVMRs))
      computedVMRs$regionType[which(computedVMRs$proportion1 < computedVMRs$proportion2)] <- "gain"
      computedVMRs$direction <- rep(-1, length(computedVMRs))
      computedVMRs$direction[which(computedVMRs$proportion1 < computedVMRs$proportion2)] <- 1
      computedVMRs <- computedVMRs[order(computedVMRs)]
    } else{
      computedVMRs <- GRanges()
    }

    if(length(computedVMRs) > 0){
      cat("Identifying VMRs...\n")
      computedVMRs <- computedVMRs[which(!is.na(computedVMRs$f_pvalue))] ## REMOVE this part if not want to discard NA value
      computedVMRs$wilcox_pvalue <- stats::p.adjust(computedVMRs$wilcox_pvalue, method = "fdr")
      computedVMRs$f_pvalue <- stats::p.adjust(computedVMRs$f_pvalue, method = "fdr")

      bufferIndex <- !is.na(computedVMRs$f_pvalue) &
        abs(computedVMRs$proportion1 - computedVMRs$proportion2) >= minProportionDifference &
        computedVMRs$sumReadsN1/computedVMRs$cytosinesCount >= minReadsPerCytosine &
        computedVMRs$sumReadsN2/computedVMRs$cytosinesCount >= minReadsPerCytosine &
        computedVMRs$cytosinesCount >= minCytosinesCount
      if (ciExcludesOne == TRUE){
        bufferIndex <- bufferIndex & (sapply(computedVMRs$F_test_result, function(x) x$conf.int[1]) > 1 |
                                        sapply(computedVMRs$F_test_result, function(x) x$conf.int[2]) < 1)
      } else {
        bufferIndex <- bufferIndex & !(sapply(computedVMRs$F_test_result, function(x) x$conf.int[1]) > 1 |
                                         sapply(computedVMRs$F_test_result, function(x) x$conf.int[2]) < 1)
      }

      if (!is.null(varRatioFc)){
        bufferIndex <- bufferIndex & (computedVMRs$var_ratio >= varRatioFc | computedVMRs$var_ratio <= 1/varRatioFc)
      }

      computedVMRs <- computedVMRs[bufferIndex]

      computedVMRs$is_DMR <- computedVMRs$wilcox_pvalue <= pValueThreshold
      computedVMRs$is_VMR <- computedVMRs$f_pvalue <= pValueThreshold
    }
  } else{
    computedVMRs <- GRanges()
  }

  return(computedVMRs)

}

#' Performs the analysis in all regions in a \code{\link[GenomicRanges]{GRanges}} object,
#' computing per-region summary statistics and statistical tests.
#'
#' @title Analyse reads inside regions (with VMR testing)
#'
#' @param methylationData A \code{GRanges} of methylation calls. Must include metadata:
#'   \describe{
#'     \item{readsM1, readsN1}{methylated and total read counts in condition 1}
#'     \item{readsM2, readsN2}{methylated and total read counts in condition 2}
#'     \item{ONT_Cm, ONT_C}{per-read methylation and coverage indices (character vectors)}
#'   }
#'
#' @param regions A \code{GRanges} of genomic intervals to summarise.
#'
#' @return A \code{GRanges} with the same ranges as \code{regions}, plus these metadata:
#' \describe{
#'   \item{sumReadsM1}{total methylated reads in condition 1}
#'   \item{sumReadsN1}{total reads in condition 1}
#'   \item{proportion1}{methylation proportion (sumReadsM1/sumReadsN1)}
#'   \item{variance1}{variance of per-read methylation proportions in condition 1}
#'   \item{sumReadsM2}{total methylated reads in condition 2}
#'   \item{sumReadsN2}{total reads in condition 2}
#'   \item{proportion2}{methylation proportion (sumReadsM2/sumReadsN2)}
#'   \item{variance2}{variance of per-read methylation proportions in condition 2}
#'   \item{cytosinesCount}{number of cytosines observed in each region}
#'   \item{wilcox_pvalue}{p-value from Wilcoxon rank-sum test comparing per-read proportions}
#'   \item{f_pvalue}{p-value from F-test comparing variances of per-read proportions}
#'   \item{var_ratio}{Ratio of variances (variance1 / variance2)}
#'   \item{wilcox_result}{Full \code{htest} object returned by \code{wilcox.test}}
#'   \item{F_test_result}{Full \code{htest} object returned by \code{var.test}}
#' }
#'
#' @author Radu Zabet and Young Jun Kim
.analyseReadsInsideRegionsVMR <- function(methylationData, regions){

  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))

  regions$sumReadsM1 <- rep(0, times=length(regions))
  regions$sumReadsN1 <- rep(0, times=length(regions))
  regions$proportion1 <- rep(0, times=length(regions))
  regions$variance1 <- rep(0, times=length(regions))
  regions$sumReadsM2 <- rep(0, times=length(regions))
  regions$sumReadsN2 <- rep(0, times=length(regions))
  regions$proportion2 <- rep(0, times=length(regions))
  regions$variance2 <- rep(0, times=length(regions))
  regions$cytosinesCount <- rep(0, times=length(regions))
  regions$wilcox_pvalue <- rep(0, times=length(regions))
  regions$f_pvalue <- rep(0, times=length(regions))
  ### run the wilcoxon test for comparing the proportion per reads between two dataset
  if (length(regionsIndexes) > 0){
    regions$sumReadsM1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM1)
    regions$sumReadsN1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN1)
    regions$sumReadsM2[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM2)
    regions$sumReadsN2[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN2)
    regions$cytosinesCount[regionsIndexes] <- sapply(methylationDataContextList,length)

    valid <- regions$cytosinesCount[regionsIndexes] > 0
    regions$proportion1[regionsIndexes[valid]] <- regions$sumReadsM1[regionsIndexes[valid]]/regions$sumReadsN1[regionsIndexes[valid]]
    regions$proportion2[regionsIndexes[valid]] <- regions$sumReadsM2[regionsIndexes[valid]]/regions$sumReadsN2[regionsIndexes[valid]]

    # get the wilcoxon test and f-test result
    test_result_list <- lapply(methylationDataContextList, .wilcox_ftestPerRead)

    wilcox_pvals <- sapply(test_result_list, function(res) res$wilcoxTest$p.value)
    var1        <- sapply(test_result_list, function(res) res$varience1)
    var2        <- sapply(test_result_list, function(res) res$varience2)
    f_pvals     <- sapply(test_result_list, function(res) {
      if (!is.null(res$fTest)) res$fTest$p.value else NA_real_
    })
    wilcox_objs <- lapply(test_result_list, function(res) res$wilcoxTest)
    f_objs     <- lapply(test_result_list, function(res) res$fTest)

    regions$wilcox_pvalue[regionsIndexes] <- wilcox_pvals
    regions$variance1[regionsIndexes]        <- var1
    regions$variance2[regionsIndexes]        <- var2
    regions$var_ratio[regionsIndexes]        <- var1 / var2
    regions$f_pvalue[regionsIndexes]         <- f_pvals
    regions$wilcox_result[regionsIndexes]    <- wilcox_objs
    regions$F_test_result[regionsIndexes]    <- f_objs
  }
  return(regions)
}


.sumReadsM1 <- function(methylationData){
  return(sum(methylationData$readsM1))
}

.sumReadsN1 <- function(methylationData){
  return(sum(methylationData$readsN1))
}

.sumReadsM2 <- function(methylationData){
  return(sum(methylationData$readsM2))
}

.sumReadsN2 <- function(methylationData){
  return(sum(methylationData$readsN2))
}

.wilcox_ftestPerRead <- function(methylationData1){
  results <- list()
  # collect the sequence index from GRanges (ONT_Cm, ONT_C)
  read_Cm_idx_list <- strsplit(unlist(methylationData1$ONT_Cm),c("_"))
  read_C_idx_list <- strsplit(unlist(methylationData1$ONT_C),c("_"))

  sample1_read_Cm_idx <- sapply(read_Cm_idx_list, function(x) x[1] == "Sample1")
  sample2_read_Cm_idx <- sapply(read_Cm_idx_list, function(x) x[1] == "Sample2")
  read_Cm_idx1_list <- read_Cm_idx_list[sample1_read_Cm_idx]
  read_Cm_idx2_list <- read_Cm_idx_list[sample2_read_Cm_idx]

  sample1_read_C_idx <- sapply(read_C_idx_list, function(x) x[1] == "Sample1")
  sample2_read_C_idx <- sapply(read_C_idx_list, function(x) x[1] == "Sample2")
  read_C_idx1_list <- read_C_idx_list[sample1_read_C_idx]
  read_C_idx2_list <- read_C_idx_list[sample2_read_C_idx]

  .make_id_vector <- function(list) {
    # list is e.g. read_Cm_idx1_list
    if (length(list)==0) return(character(0))
    # extract the "second field", coerce to character, drop any NAs
    tmp <- as.character(sapply(list, function(x) x[2], USE.NAMES=FALSE))
    tmp <- tmp[!is.na(tmp)]
    if (length(tmp)==0) return(character(0))
    unique(unlist(strsplit(tmp, ",", fixed=TRUE), use.names=FALSE))
  }
  read_Cm_idx1 <- .make_id_vector(read_Cm_idx1_list)
  read_Cm_idx2 <- .make_id_vector(read_Cm_idx2_list)

  read_C_idx1  <- .make_id_vector(read_C_idx1_list)
  read_C_idx2  <- .make_id_vector(read_C_idx2_list)

  total_idx1 <- sort(unique(c(read_Cm_idx1,read_C_idx1)))
  total_idx2 <- sort(unique(c(read_Cm_idx2,read_C_idx2)))

  # compute proportions per reads using the GRanges (ONT_Cm, ONT_C)
  proportions_S1 <- rep(0, length(total_idx1))
  proportions_S2 <- rep(0, length(total_idx1))

  for (k in seq_along(total_idx1)){
    id <- total_idx1[k]
    per_read_Cm_S1 <- sapply(mcols(methylationData1)$ONT_Cm, function(x) {
      s1 <- x[grepl("^Sample1_", x)]
      if (length(s1) == 0) return(0L)
      sum(str_count(s1, fixed(id)))
    })

    per_read_C_S1 <- sapply(mcols(methylationData1)$ONT_C, function(x) {
      s1 <- x[grepl("^Sample1_", x)]
      if (length(s1) == 0) return(0L)
      sum(str_count(s1, fixed(id)))
    })
    per_read_S1_M <-sum(per_read_Cm_S1)
    per_read_S1_N <-sum(per_read_Cm_S1)+sum(per_read_C_S1)
    proportions_S1[k] <- per_read_S1_M/per_read_S1_N
  }

  for (k in seq_along(total_idx2)){
    id <- total_idx2[k]
    per_read_Cm_S2 <- sapply(mcols(methylationData1)$ONT_Cm, function(x) {
      s2 <- x[grepl("^Sample2_", x)]
      if (length(s2) == 0) return(0L)
      sum(str_count(s2, fixed(id)))
    })

    per_read_C_S2 <- sapply(mcols(methylationData1)$ONT_C, function(x) {
      s2 <- x[grepl("^Sample2_", x)]
      if (length(s2) == 0) return(0L)
      sum(str_count(s2, fixed(id)))
    })

    per_read_S2_M <-sum(per_read_Cm_S2)
    per_read_S2_N <-sum(per_read_Cm_S2)+sum(per_read_C_S2)
    proportions_S2[k] <- per_read_S2_M/per_read_S2_N
    # ## DEBUG
    # cat(c("per_read_S2_M:",per_read_S2_M," per_read_S2_N:",per_read_S2_N," per_read_S2_P:",per_read_S2_M/per_read_S2_N,"\n"))
  }

  # run Wilcoxon test and f-test
  if (length(proportions_S1) > 0 & length(proportions_S2) > 0) {
    wt <- try(wilcox.test(proportions_S1, proportions_S2), silent=TRUE)
    if (inherits(wt, "try-error") || !inherits(wt, "htest")) {
      wt <- structure(
        list(
          statistic   = setNames(NA_real_, "W"),
          p.value     = NA_real_,
          alternative = "true location shift is not equal to 0",
          method      = "Wilcoxon rank sum test",
          data.name   = paste0("proportions_S1 and proportions_S2")
        ),
        class = "htest"
      )
    }
  }
  results$wilcoxTest<- wt
  results$varience1 <- var(proportions_S1)
  results$varience2 <- var(proportions_S2)

  if(length(proportions_S1)>=2 & length(proportions_S2)>=2 &
     var(proportions_S1, na.rm=TRUE) > 0 & var(proportions_S2, na.rm=TRUE) > 0){
    results$fTest <- var.test(proportions_S1,proportions_S2)
  } else {
    results$fTest <- NULL
  }
 return(results)
}

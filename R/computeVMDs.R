#' This function computes the variance methylated domains between pre-set
#' min and max proportion values.  
#'
#' @title Compute VMDs
#' @param methylationData the methylation data in condition
#' (see \code{\link{ontSampleGRangesList}}).
#' @param regions a \code{\link{GRanges}} object with the regions where to 
#' compute the VMDs. If \code{NULL}, the VMDs are computed genome-wide.
#' @param context the context in which the VMDs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param binSize the size of the tiling bins in nucleotides. This parameter is 
#' required only if the selected method is \code{"bins"}.
#' @param minCytosinesCount VMDs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param sdCutoffMethod Character string specifying how to determine the cutoff
#' for filtering VMDs based on their methylation variance (standard deviation).
#' Available options are:
#' \describe{
#'   \item{\code{"per.high"}}{Selects the top \code{percentage} of regions with the highest variance (standard deviation).}
#'   \item{\code{"per.low"}}{Selects the bottom \code{percentage} of regions with the lowest variance.}
#'   \item{\code{"EDE.high"}}{Uses the elbow point (inflection/knee) from the descendingly sorted variance values to determine a data-driven high-variance cutoff. Retains regions with SD above this elbow point.}
#'   \item{\code{"EDE.low"}}{Uses the elbow point from the ascendingly sorted variance values to define a low-variance cutoff. Retains regions with SD below this point.}
#' }
#' This allows either quantile-based filtering or automatic detection of variance thresholds based on distribution shape.
#' @param percentage Numeric cutoff used when \code{sdCutoffMethod} is set to
#' \code{"per.high"} or \code{"per.low"}. Represents the quantile threshold:
#' for example, \code{percentage = 0.05} keeps the top 5\% or bottom 5\% of
#' bins based on standard deviation, depending on the selected method.
#' @param minGap VMDs separated by a gap of at least \code{minGap} are not 
#' merged. Note that only VMDs where the change in methylation is in the same 
#' direction are joined.
#' @param minSize VMDs with a size smaller than \code{minSize} are discarded.
#' @param minReadsPerCytosine  VMDs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return the VMDs stored as a \code{\link{GRanges}} object with the following 
#' metadata columns:
#' \describe{
#'  \item{context}{the context in which the VMDs was computed (\code{"CG"}, 
#'  \code{"CHG"} or \code{"CHH"}).}
#'  \item{sumReadsM}{the number of methylated reads.}
#'  \item{sumReadsN}{the total number of reads.} 
#'  \item{proportion}{the proportion from total methylated reads.} 
#'  \item{cytosinesCount}{the number of cytosines in the VMDs.}
#'  \item{mean}{mean value comparing per‐read proportions}
#'  \item{sd}{standard deviation comparing per‐read proportions}
#'  \item{w_mean}{weighted mean value comparing per‐read proportions}
#'  \item{w_sd}{weighted standard deviation comparing per‐read proportions} 
#' }
#' @seealso \code{\link{readONTbam}}, \code{\link{filterVMDs}}
#' and \code{\link{analyseReadsInsideRegionsForCondition}}
#' @examples
#' \dontrun{
#' # load the ONT methylation data 
#' data(ontSampleGRangesList)
#' 
#' # the regions where to compute the VMDs
#' chr1_ranges <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,1E6+6E5))
#' 
#' # compute the VMDs in CG context with bins method
#' VMDsBinsCG <- computeVMDs(ontSampleGRangesList[["GM18501"]],
#'                           regions = NULL,
#'                           context = "CG",
#'                           binSize = 100,
#'                           minCytosinesCount = 4,
#'                           sdCutoffMethod = "EDE.high",
#'                           percentage = 0.05,
#'                           minGap = 200,
#'                           minSize = 50,
#'                           minReadsPerCytosine = 4,
#'                           parallel = FALSE,
#'                           BPPARAM = NULL,
#'                           cores = 1)
#' }
#' @author Nicolae Radu Zabet, Jonathan Michael Foonlan Tsang and Young Jun Kim
#' @import inflection 
#' @export
 

computeVMDs <- function(methylationData, 
                        regions = NULL, 
                        context = "CG", 
                        binSize = 100,
                        minCytosinesCount = 4, 
                        sdCutoffMethod = "per.high",
                        percentage = 0.05,
                        minGap = 200, 
                        minSize = 50, 
                        minReadsPerCytosine = 4, 
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
  if (!is.null(cores)){
    .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                  " the number of cores to used when computing the DMRs needs to be an integer  hirger or equal to 1.")
    .stopIfNotAll(BPPARAM$workers >= cores, paste0("the cores should be smaller than system's cores, current system cores: ",BPPARAM$workers))
    BPPARAM$workers <- cores
  } else {
    cores <- BPPARAM$workers
  }
  cat("Current parallel setting, BPPARAM: ", capture.output(BPPARAM),sep = "\n")
  
  .validateMethylationData(methylationData, variableName="methylationData")
  
  regions <- .validateGRanges(regions, methylationData)
  
  .validateContext(context)
  
  .stopIfNotAll(c(.isInteger(binSize, positive=TRUE)),
                " the bin size used by the method is an integer higher than 0")
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(sdCutoffMethod),
                  all(is.character(sdCutoffMethod)),
                  length(sdCutoffMethod) == 1,
                  all(sdCutoffMethod %in% c("per.high","per.low","EDE.high","EDE.low"))),
                " sdCutoffMethod can be only per.high, per.low, EDE.high or EDE.low")
  
  if (sdCutoffMethod %in% c("per.high", "per.low")){
    .stopIfNotAll(c(!is.null(percentage),is.numeric(percentage), percentage >= 0, percentage <= 1),
                  " percentage needs to be a numeric value between 0 and 1.")
  }
  
  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between VMDs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minSize, positive=TRUE)),
                " the minimum size of a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")
  
  computedVMDs <- GRanges()
  
  computedVMDs <- .computeVMDsBins(methylationData = methylationData, 
                                   regions = regions, 
                                   context = context, 
                                   binSize = binSize,
                                   minCytosinesCount = minCytosinesCount, 
                                   minGap = minGap, 
                                   minSize = minSize, 
                                   minReadsPerCytosine = minReadsPerCytosine,
                                   sdCutoffMethod = sdCutoffMethod,
                                   percentage = percentage,
                                   cores = cores,
                                   parallel = parallel,
                                   BPPARAM = BPPARAM)
  
  return(computedVMDs)  
}


#' This function computes the variance methylated regions 
#' using the bins method. 
.computeVMDsBins <- function(methylationData,
                             regions = NULL, 
                             context = "CG", 
                             binSize = 100,
                             minCytosinesCount = 4,
                             minGap = 200, 
                             minSize = 50, 
                             minReadsPerCytosine = 4, 
                             sdCutoffMethod = "per.high",
                             percentage = 0.05,
                             cores = 1,
                             parallel = parallel,
                             BPPARAM = BPPARAM) {
  
  regions <- reduce(regions)
  
  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  
  contextMethylationData <- methylationData[methylationData$context%in%context]
  rm(methylationData)
  localContextMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, regions))]
  rm(contextMethylationData)  
  
  regionsList <- .splitGRangesEqualy(regions, cores)
  
  # inner loop function for BiocParallel::bplapply
  .computeVMDsBinsLoop = function(i){
    computedVMDs <- GRanges()
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]
      
      
      cat("Computing VMDs at ",.printGenomicRanges(currentRegion),"\n")
      
      seqs <- seq(start(currentRegion), (end(currentRegion)-binSize), by = binSize);
      
      bins <- GRanges(seqnames(currentRegion), IRanges(seqs, (seqs+binSize-1)))
      
      overlapsBins <- findOverlaps(localContextMethylationData, currentRegion)
      
      if(length(overlapsBins) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlapsBins)]
        
        cat("Count inside each bin...\n")
        #bins <- .analyseReadsInsideRegions(localMethylationData, bins, context, cores)
        bins <- .analyseReadsInsideBinsVMDs(localMethylationData, bins, currentRegion)
        
        cat("Filter the bins...\n")
        # Get rid of the bins with fewer than minCytosinesCount cytosines inside.  
        bins  <- bins[bins$cytosinesCount >= minCytosinesCount]
        
        # Get rid of the bins with fewer than minReadsPerCytosine reads per cytosine.  
        bins  <- bins[(bins$sumReadsN/bins$cytosinesCount >= minReadsPerCytosine)]
        
        # filtering VMD in this part
        
        cat("Identifying VMDs...\n")    
        bins$context <- rep(paste(context, collapse = "_"), length(bins))
        bins$direction <- rep(NA, length(bins))
        bins$direction <- sign(bins$proportion)
        
        # Select the crude list of VMDs
        VMDs <- bins[!is.na(bins$direction) & bins$direction == 1]
        
        # append current VMDs to the global list of VMDs 
        if(length(computedVMDs) == 0){
          computedVMDs <- VMDs
        } else{
          computedVMDs <- c(computedVMDs,VMDs)
        }
      }
    }
    return(computedVMDs)
  }
  
  # compute the VMDs
  if(cores > 1){
    cat("Compute the VMDs using ", cores, "cores\n")
    computedVMDs <- BiocParallel::bplapply(1:length(regionsList), .computeVMDsBinsLoop, BPPARAM = BPPARAM)
  } else {
    computedVMDs <- lapply(1:length(regionsList), .computeVMDsBinsLoop)
  }
  
  computedVMDs <- unlist(GRangesList(computedVMDs))
  
  if(length(computedVMDs) > 0){
    
    cat("Merge adjacent VMDs\n")    
    computedVMDs <- computedVMDs[order(computedVMDs)]
    
    cat("Merge VMDs iteratively\n")    
    # Get rid of small gaps between VMDs.
    if(minGap > 0){
      computedVMDs <- .smartMergePMDs(computedVMDs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minReadsPerCytosine = minReadsPerCytosine,
                                      minMethylation = 0,
                                      maxMethylation = 1,
                                      cores = cores,
                                      BPPARAM = BPPARAM)
    }  
    computedVMDs <- computedVMDs[!duplicated(granges(computedVMDs))]
    hits <- findOverlaps(computedVMDs, computedVMDs, type = "within", ignore.strand = FALSE)
    computedVMDs <- computedVMDs[-unique(queryHits(hits)[queryHits(hits) != subjectHits(hits)])]
    
    cat("Filter VMDs \n") 
    
    computedVMDs <- DMRcaller::filterVMDs(localContextMethylationData,
                               computedVMDs,
                               context = context,
                               minCytosinesCount = minCytosinesCount,
                               minReadsPerCytosine = minReadsPerCytosine,
                               sdCutoffMethod = sdCutoffMethod,
                               percentage = percentage,
                               parallel = parallel,
                               BPPARAM = BPPARAM,
                               cores = cores)
    
   
    if(length(computedVMDs) > 0){
      #remove small VMDs 
      computedVMDs <- computedVMDs[width(computedVMDs) >= minSize]
      if(length(computedVMDs) > 0){
        #remove VMDswith few cytosines
        computedVMDs <- computedVMDs[!is.na(computedVMDs$cytosinesCount) & computedVMDs$cytosinesCount >= minCytosinesCount]
      }
    }
  }  
  if (length(computedVMDs) >0 & length(computedVMDs$direction) != 0){
    computedVMDs$direction <- NULL
  }
  return(computedVMDs)
}

.analyseReadsInsideBinsVMDs <- function(methylationData, bins, currentRegion){
  
  binSize <- min(unique(width(bins)))
  #Rcpp
  readsM <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsM, windowSize = binSize)
  sumReadsM <- readsM[seq(1,length(readsM)-binSize, by=binSize)]
  
  readsN <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), methylationData$readsN, windowSize = binSize)
  sumReadsN <- readsN[seq(1,length(readsN)-binSize, by=binSize)]
  
  proportion <- sumReadsM/sumReadsN
  
  cytosines <- .movingSum(start(currentRegion), end(currentRegion), start(methylationData), rep(1, length(start(methylationData))), windowSize = binSize)
  cytosinesCount <- cytosines[seq(1,length(cytosines)-binSize, by=binSize)]
  
  
  overlaps <- findOverlaps(methylationData, bins, ignore.strand = TRUE)
  methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  bin_result_list <- lapply(methylationDataContextList, .mean_sdPerRead)
  
  mean1 <- sapply(bin_result_list, function(res) res$mean)
  sd1   <- sapply(bin_result_list, function(res) res$sd)
  w_mean1 <- sapply(bin_result_list, function(res) res$mean_w)
  w_sd1   <- sapply(bin_result_list, function(res) res$sd_w)
  
  bins$sumReadsM <- sumReadsM
  bins$sumReadsN <- sumReadsN   
  bins$proportion <- proportion        
  bins$cytosinesCount <- cytosinesCount 
  
  bins$mean[regionsIndexes] <- mean1
  bins$sd[regionsIndexes]   <- sd1
  bins$w_mean[regionsIndexes] <- w_mean1
  bins$w_sd[regionsIndexes]   <- w_sd1
  
  return(bins)
}


#' This function verifies whether a set of potential VMDs (e.g. genes, 
#' transposons, CpG islands) are variance methylated or not.
#'
#' @title Filter VMDs 
#' @param methylationData the methylation data in condition
#' (see \code{\link{ontSampleGRangesList}}).
#' @param potentialVMDs a \code{\link{GRanges}} object with potential VMDs 
#' where to compute the VMDs. This can be a a list of gene and/or transposable 
#' elements coordinates.
#' @param context the context in which the VMDs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param minCytosinesCount VMDs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param minReadsPerCytosine  VMDs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param sdCutoffMethod Character string specifying how to determine the cutoff
#' for filtering VMDs based on their methylation variance (weighted standard deviation).
#' Available options are:
#' \describe{
#'   \item{\code{"per.high"}}{Selects the top \code{percentage} of regions with the highest variance (standard deviation).}
#'   \item{\code{"per.low"}}{Selects the bottom \code{percentage} of regions with the lowest variance.}
#'   \item{\code{"EDE.high"}}{Uses the elbow point (inflection/knee) from the descendingly sorted variance values to determine a data-driven high-variance cutoff. Retains regions with SD above this elbow point.}
#'   \item{\code{"EDE.low"}}{Uses the elbow point from the ascendingly sorted variance values to define a low-variance cutoff. Retains regions with SD below this point.}
#' }
#' This allows either quantile-based filtering or automatic detection of variance thresholds based on distribution shape.
#' @param percentage Numeric cutoff used when \code{sdCutoffMethod} is set to
#' \code{"per.high"} or \code{"per.low"}. Represents the quantile threshold:
#' for example, \code{percentage = 0.05} keeps the top 5\% or bottom 5\% of
#' bins based on weighted standard deviation, depending on the selected method.
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return a \code{\link{GRanges}} object with 9 metadata columns that contain 
#' the VMDs; see \code{\link{computeVMDs}}.
#' @seealso \code{\link{computeVMDs}} 
#' and \code{\link{analyseReadsInsideRegionsForCondition}}  
#' @examples
#' # load the ONT methylation data 
#' data(ontSampleGRangesList)
#' # load the gene annotation data
#' data(GEs_hg38)
#' 
#' # select the transcript
#' transcript <- GEs_hg38[which(GEs_hg38$type == "transcript")]
#' 
#' # the regions where to compute the VMDs
#' regions <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' transcript <- transcript[overlapsAny(transcript, regions)]
#' 
#' # filter genes that are variance methylated in the two conditions
#' VMDsGenesCG <- filterVMDs(ontSampleGRangesList[["GM18501"]], 
#'                potentialVMDs = transcript, 
#'                context = "CG", sdCutoffMethod = "per.high", percentage = 0.05,
#'                minCytosinesCount = 4, minReadsPerCytosine = 3, cores = 1)
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @export
filterVMDs <- function(methylationData, 
                       potentialVMDs, 
                       context = "CG",   
                       minCytosinesCount = 4, 
                       minReadsPerCytosine = 3, 
                       sdCutoffMethod = "per.high",
                       percentage = 0.05,
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
  if (!is.null(cores)){
    .stopIfNotAll(c(.isInteger(cores, positive=TRUE)), 
                  " the number of cores to used when computing the DMRs needs to be an integer  hirger or equal to 1.")
    .stopIfNotAll(BPPARAM$workers >= cores, paste0("the cores should be smaller than system's cores, current system cores: ",BPPARAM$workers))
    BPPARAM$workers <- cores
  } else {
    cores <- BPPARAM$workers
  }
  cat("Current parallel setting, BPPARAM: ", capture.output(BPPARAM),sep = "\n")
  
  .validateMethylationData(methylationData, variableName="methylationData")
  
  regions <- getWholeChromosomes(methylationData)
  
  .validateContext(context)
  
  .validateGRanges(potentialVMDs, generateGenomeWide=FALSE, variableName="potentialVMDs", minLength=NULL)
  
  regions <- .validateGRanges(regions, methylationData)
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum number of reads in a bin is an integer higher or equal to 0")
  
  regions <- reduce(regions)
  
  if(length(potentialVMDs) > 0){
    
    
    
    # extract the methylation data in the correct context
    cat("Extract methylation in the corresponding context \n")
    
    contextMethylationData <- methylationData[methylationData$context%in%context]
    rm(methylationData)
    localContextMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, regions))]
    rm(contextMethylationData)  
    
    regionsList <- .splitGRangesEqualy(regions, cores)
    
    # inner loop function for BiocParallel::bplapply
    .filterVMDsLoop = function(i){
      computedVMDs <- GRanges()  
      for(index in 1:length(regionsList[[i]])){
        currentRegion <- regionsList[[i]][index]
        
        
        cat("Computing VMDs at ",.printGenomicRanges(currentRegion),"\n")
        
        cat("Selecting data...\n")
        
        # Select the points in methylationData that we're interested in. These are the 
        # points that lie within 'regions', as well as any that lie within 
        # window.size of them. 
        
        overlapsPotentialVMDs <- findOverlaps(potentialVMDs, currentRegion)
        if(length(overlapsPotentialVMDs) > 0){
          potentialVMDsLocal <- potentialVMDs[queryHits(overlapsPotentialVMDs)]
          
          localMethylationData <- localContextMethylationData[queryHits(findOverlaps(localContextMethylationData, currentRegion))]
          potentialVMDsLocal <- .analyseReadsInsideRegionsVMDs(localMethylationData, potentialVMDsLocal)
          
          if(length(computedVMDs) == 0){
            computedVMDs <- potentialVMDsLocal
          } else{
            computedVMDs <- c(computedVMDs,potentialVMDsLocal)
          }
        } 
      }
      return(computedVMDs)
    }
    
    # compute the VMDs
    if(cores > 1){
      cat("Compute the VMDs using ", cores, "cores\n")
      computedVMDs <- BiocParallel::bplapply(1:length(regionsList), .filterVMDsLoop, BPPARAM = BPPARAM)
    } else {
      computedVMDs <- lapply(1:length(regionsList), .filterVMDsLoop)
    }
    
    
    computedVMDs <-  unlist(GRangesList(computedVMDs))
    
    if(length(computedVMDs) > 0){
      cat("Identifying VMDs...\n")    
      
      bufferIndex <- !is.na(computedVMDs$proportion) &
        computedVMDs$sumReadsN/computedVMDs$cytosinesCount >= minReadsPerCytosine &
        computedVMDs$cytosinesCount >= minCytosinesCount
      
      computedVMDs <- computedVMDs[bufferIndex]  
    }
  } else{
    computedVMDs <- GRanges() 
  }
  
  if(length(computedVMDs) > 0){
    computedVMDs <- computedVMDs[order(computedVMDs)]
  }
  ## cutoff by percentage method
  if (sdCutoffMethod %in% c("per.high","per.low")){
    # # Convert VMD result to data frame
    # df <- as.data.frame(computedVMDs)
    
    # # Sort by SD in descending order
    # df_sorted <- df[order(-df$w_sd), ]
    # 
    # # Create x-axis index
    # df_sorted$bin_index <- seq_len(nrow(df_sorted))
    
    # Setting the cutoff of sd
    cutoff_high <- quantile(computedVMDs$w_sd, 1-percentage, na.rm = TRUE)
    cutoff_low <- quantile(computedVMDs$w_sd, percentage, na.rm = TRUE)
    
    # col <- c("#D55E00","#E69F00", "#0072B2", "#56B4E9", "#F0E442", "#009E73")
    # # Plot for visual inspection 
    # plot(df_sorted$bin_index, df_sorted$w_sd, type="l", xlab="Bins (sorted by SD)", col = "black",
    #      ylab="Standard Deviation of Methylation", main = "SD distribution across bins")
    # abline(h = cutoff_high, lty = 3, col = col[1])  # horizontal dashed line at SD threshold
    # abline(h = cutoff_low, lty = 3, col = col[6])
    # legend(x = "topright", legend=c("cutoff_high", "cutoff_low"), lty = c(3,3),
    #        col = c(col[1],col[6]),   bty = "n",  cex = 0.7,inset = c(0.01, 0.01)
    # )
    if (sdCutoffMethod == "per.high"){
      computedVMDs <- computedVMDs[computedVMDs$w_sd >= cutoff_high]
    } else if (sdCutoffMethod == "per.low"){
      computedVMDs <- computedVMDs[computedVMDs$w_sd <= cutoff_low]
    }
  }
  
  ## cutoff by EDE method
  if (sdCutoffMethod %in% c("EDE.high","EDE.low")){
    # Sort them in decreasing order so the largest variance bins come first
    y <- sort(computedVMDs$w_sd, decreasing = TRUE)   
    x <- seq_along(y)   
    
    # Find the elbow point (knee) in the curve
    ede_res <- inflection::ede(x, y, 0)
    k_idx2 <- ede_res[1]
    k_idx3 <- ede_res[2]
    k_idx4 <- ede_res[3]
    knee_y2 <- y[k_idx2]
    knee_y3 <- y[k_idx3]
    knee_y4 <- y[k_idx4]
    
    # col <- c("#D55E00","#E69F00", "#0072B2", "#56B4E9", "#F0E442", "#009E73")
    # # Plot for visual inspection 
    # plot(x, y, type="l", xlab="Bins (sorted by SD)", col = "black",
    #      ylab="Standard Deviation of Methylation", main = "SD distribution across bins")
    # abline(h = knee_y2, lty = 3, col = col[2])  # horizontal dashed line at SD threshold
    # abline(h = knee_y4, lty = 3, col = col[3])
    # abline(h = knee_y3, lty = 3, col = col[4])
    # legend(x = "topright", legend=c("EDE j1", "EDE j2", "EDE j3"), lty = c(3,3,3),
    #        col = col[2:4],   bty = "n",  cex = 0.7,inset = c(0.01, 0.01)
    # )
    if (sdCutoffMethod == "EDE.high"){
      computedVMDs <- computedVMDs[computedVMDs$w_sd >= knee_y2]
    } else if (sdCutoffMethod == "EDE.low"){
      computedVMDs <- computedVMDs[computedVMDs$w_sd <= knee_y3]
    }
  }
  return(computedVMDs)
  
}

#' Performs the analysis in all regions in a \code{\link{GRanges}} object,
#' computing per-region mean and standard deviation of proportions.
#'
#' @title Analyse reads inside regions with computing mean and standard deviation per-read proportions
#'
#' @param methylationData A \code{GRanges} of methylation calls. Must include metadata:
#'   \describe{
#'     \item{readsM, readsN}{methylated and total read counts in condition}
#'     \item{ONT_Cm, ONT_C}{per-read methylation and coverage indices (character vectors)}
#'   }
#'
#' @param regions A \code{GRanges} of genomic intervals to summarise.
#'
#' @return A \code{GRanges} with the same ranges as \code{regions}, plus these metadata:
#' \describe{
#'   \item{sumReadsM}{total methylated reads in condition}
#'   \item{sumReadsN}{total reads in condition}
#'   \item{proportion}{methylation proportion (sumReadsM/sumReadsN)}
#'   \item{cytosinesCount}{number of cytosines observed in each region}
#'   \item{mean}{mean value comparing per‐read proportions}
#'   \item{sd}{standard deviation comparing per‐read proportions}
#'   \item{w_mean}{weighted mean value comparing per‐read proportions}
#'   \item{w_sd}{weighted standard deviation comparing per‐read proportions}
#' } 
#'       
#' @author Radu Zabet and Young Jun Kim
.analyseReadsInsideRegionsVMDs <- function(methylationData, regions){
  
  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- S4Vectors::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))
  
  regions$sumReadsM <- rep(0, times=length(regions))
  regions$sumReadsN <- rep(0, times=length(regions))    
  regions$proportion <- rep(0, times=length(regions)) 
  regions$cytosinesCount <- rep(0, times=length(regions))
  regions$mean <- rep(0, times=length(regions))
  regions$sd <- rep(0, times=length(regions))
  regions$w_mean <- rep(0, times=length(regions))
  regions$w_sd <- rep(0, times=length(regions))
  
  
  if(length(regionsIndexes) > 0){  
    regions$sumReadsM[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsM)
    regions$sumReadsN[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsN)               
    regions$cytosinesCount[regionsIndexes] <- sapply(methylationDataContextList,length)
    
    valid <- regions$cytosinesCount[regionsIndexes] > 0
    regions$proportion[regionsIndexes[valid]] <- regions$sumReadsM[regionsIndexes[valid]]/regions$sumReadsN[regionsIndexes[valid]]
    
    # get the mean and sd result per reads
    test_result_list <- lapply(methylationDataContextList, .mean_sdPerRead)
    
    mean1 <- sapply(test_result_list, function(res) res$mean)
    sd1   <- sapply(test_result_list, function(res) res$sd)
    w_mean1 <- sapply(test_result_list, function(res) res$mean_w)
    w_sd1   <- sapply(test_result_list, function(res) res$sd_w)
    
    regions$mean[regionsIndexes] <- mean1
    regions$sd[regionsIndexes]   <- sd1
    regions$w_mean[regionsIndexes] <- w_mean1
    regions$w_sd[regionsIndexes]   <- w_sd1
  }
  return(regions)
}


.sumReadsM <- function(methylationData){
  return(sum(methylationData$readsM))
}
.sumReadsN <- function(methylationData){
  return(sum(methylationData$readsN))
}

.mean_sdPerRead <- function(methylationData1){
  if (!("ONT_Cm" %in% names(mcols(methylationData1))) |
      !("ONT_C"  %in% names(mcols(methylationData1)))) {
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
  
  read_Cm_ids <- split_ids(unlist(methylationData1$ONT_Cm))
  read_C_ids  <- split_ids(unlist(methylationData1$ONT_C))
  total_idx   <- sort(unique(c(read_Cm_ids, read_C_ids)))
  total_idx   <- total_idx[nzchar(total_idx)]
  
  if (length(total_idx) == 0L) {
    return(list(mean = NA_real_, sd = NA_real_, mean_w = NA_real_, sd_w = NA_real_))
  }

  Cm_list <- lapply(mcols(methylationData1)$ONT_Cm, split_ids)
  C_list  <- lapply(mcols(methylationData1)$ONT_C,  split_ids)
  
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
  p  <- proportions_S[ok]
  w  <- N_vec[ok]
  
  
  if (length(p) == 0L) {
    return(list(mean = NA_real_, sd = NA_real_, mean_w = NA_real_, sd_w = NA_real_))
  }
  
  mean_u <- mean(p)
  sd_u   <- sd(p)
  
  if (length(w) > 0 && sum(w) > 0) {
    mu_w  <- sum(w * p) / sum(w)
    var_w <- sum(w * (p - mu_w)^2) / sum(w)   
    sd_w  <- sqrt(var_w)
  } else {
    mu_w <- NA_real_
    sd_w <- NA_real_
  }
  
  list(mean = mean_u, sd = sd_u, mean_w = mu_w, sd_w = sd_w)
}

#' This function computes the partially methylated domains between pre-set
#' min and max proportion values.  
#'
#' @title Compute PMDs
#' @param methylationData the methylation data in condition
#' (see \code{\link{ontSampleGRangesList}}).
#' @param regions a \code{\link[GenomicRanges]{GRanges}} object with the regions where to 
#' compute the PMDs. If \code{NULL}, the PMDs are computed genome-wide.
#' @param context the context in which the PMDs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param method Character string specifying the algorithm for PMD detection. 
#' If \code{"noise_filter"}, a sliding window of size \code{windowSize} is applied 
#' with the specified \code{kernelFunction} (and \code{lambda} for a Gaussian kernel) 
#' to smooth methylation counts before calling and merging PMDs. 
#' If \code{"neighbourhood"}, individual partially methylated cytosines are 
#' identified first and then merged into PMDs. 
#' If \code{"bins"}, the genome is partitioned into fixed bins of size \code{binSize}, 
#' partially methylation is sorted per bin, and significant bins are merged. 
#' @param windowSize the size of the triangle base measured in nucleotides. 
#' This parameter is required only if the selected method is 
#' \code{"noise_filter"}. 
#' @param kernelFunction a \code{character} indicating which kernel function to 
#' be used. Can be one of \code{"uniform"}, \code{"triangular"}, 
#' \code{"gaussian"} or \code{"epanechnicov"}. This is required only if the 
#' selected method is \code{"noise_filter"}. 
#' @param lambda numeric value required for the Gaussian filter 
#' (\code{K(x) = exp(-lambda*x^2)}). This is required only if the selected 
#' method is \code{"noise_filter"} and the selected kernel function is 
#' \code{"gaussian"}. 
#' @param binSize the size of the tiling bins in nucleotides. This parameter is 
#' required only if the selected method is \code{"bins"}.
#' @param minCytosinesCount PMDs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param minMethylation Numeric [0,1]; minimum mean methylation proportion.
#' @param maxMethylation Numeric [0,1]; maximum mean methylation proportion.
#' @param minGap PMDs separated by a gap of at least \code{minGap} are not 
#' merged. Note that only PMDs where the change in methylation is in the same 
#' direction are joined.
#' @param minSize PMDs with a size smaller than \code{minSize} are discarded.
#' @param minReadsPerCytosine  PMDs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return the PMDs stored as a \code{\link[GenomicRanges]{GRanges}} object with the following 
#' metadata columns:
#' \describe{
#'  \item{context}{the context in which the PMDs was computed (\code{"CG"}, 
#'  \code{"CHG"} or \code{"CHH"}).}
#'  \item{sumReadsM}{the number of methylated reads.}
#'  \item{sumReadsN}{the total number of reads.} 
#'  \item{proportion}{the proportion methylated reads filtered between
#'  \code{minMethylation} and \code{maxMethylation}.} 
#'  \item{cytosinesCount}{the number of cytosines in the PMDs.} 
#' }
#' @seealso \code{\link{readONTbam}}, \code{\link{filterPMDs}}, \code{\link{mergePMDsIteratively}}, 
#' \code{\link{analyseReadsInsideRegionsForConditionPMD}} and 
#' \code{\link{PMDsNoiseFilterCG}}
#' @examples
#' # load the ONT methylation data 
#' data(ontSampleGRangesList)
#' 
#' # the regions where to compute the PMDs
#' chr1_ranges <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' 
#' # compute the PMDs in CG context with noise_filter method
#' PMDsNoiseFilterCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                                  regions = chr1_ranges,
#'                                  context = "CG",
#'                                  windowSize = 100,
#'                                  method = "noise_filter",
#'                                  kernelFunction = "triangular", 
#'                                  lambda = 0.5,
#'                                  minCytosinesCount = 4, 
#'                                  minMethylation = 0.4, 
#'                                  maxMethylation = 0.6, 
#'                                  minGap = 200, 
#'                                  minSize = 50, 
#'                                  minReadsPerCytosine = 4, 
#'                                  cores = 1,
#'                                  parallel = FALSE)
#'  
#' \dontrun{                                
#' # compute the PMDs in CG context with neighbourhood method
#' PMDsNeighbourhoodCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                                    regions = chr1_ranges,
#'                                    context = "CG",
#'                                    method = "neighbourhood"
#'                                    minCytosinesCount = 4, 
#'                                    minMethylation = 0.4, 
#'                                    maxMethylation = 0.6, 
#'                                    minGap = 200, 
#'                                    minSize = 50, 
#'                                    minReadsPerCytosine = 4, 
#'                                    cores = 1,
#'                                    parallel = FALSE)
#'                                    
#' # compute the PMDs in CG context with bins method
#' PMDsBinsCG <- computePMDs(ontSampleGRangesList[["GM18501"]],
#'                           regions = chr1_ranges,
#'                           context = "CG",
#'                           method = "bins",
#'                           binSize = 100,
#'                           minCytosinesCount = 4, 
#'                           minMethylation = 0.4, 
#'                           maxMethylation = 0.6, 
#'                           minGap = 200, 
#'                           minSize = 50, 
#'                           minReadsPerCytosine = 4, 
#'                           cores = 1,
#'                           parallel = FALSE)
#' }
#' @author Nicolae Radu Zabet, Jonathan Michael Foonlan Tsang and Young Jun Kim
#' @export
computePMDs <- function(methylationData, 
                        regions = NULL, 
                        context = "CG", 
                        method="noise_filter",
                        windowSize = 100,
                        kernelFunction = "triangular", 
                        lambda = 0.5,
                        binSize = 100,
                        minCytosinesCount = 4, 
                        minMethylation = 0.4, 
                        maxMethylation = 0.6, 
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
  
  .validateContext(context)
  
  .stopIfNotAll(c(!is.null(method), 
                  all(is.character(method)),
                  length(method) == 1,
                  all(method %in% c("noise_filter","neighbourhood","bins"))),
                " method can be only noise_filter, neighbourhood or bins")  
  
  if(method == "noise_filter"){
    .stopIfNotAll(c(.isInteger(windowSize, positive=TRUE)), 
                  " the window size used by the interpolation method is an integer higher than 0")
    
    .stopIfNotAll(c(!is.null(kernelFunction), 
                    kernelFunction%in%c("uniform", "triangular", "gaussian", "epanechnicov")), 
                  paste("Unknown kernel function: ", kernelFunction, ". 
                        kernelFunction should be one of \"uniform\", \"triangular\", 
                        \"gaussian\", \"epanechnicov\"",sep=""))
    
    if(kernelFunction == "gaussian"){
      .stopIfNotAll(c(!is.null(lambda),
                      is.numeric(lambda)), 
                    " lambda needs to be a numeric value.")
    }
    
  }
  
  if(method == "bins"){
    .stopIfNotAll(c(.isInteger(binSize, positive=TRUE)), 
                  " the bin size used by the method is an integer higher than 0")
    
  }
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(minMethylation), is.numeric(minMethylation), minMethylation >= 0), 
                " minMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(maxMethylation), is.numeric(maxMethylation), maxMethylation >= 0), 
                " maxMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(maxMethylation >= minMethylation), 
                " maxMethylation should be higher than minMethylation value")
  
  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between PMDs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minSize, positive=TRUE)),
                " the minimum size of a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")
  
  computedPMDs <- GRanges()
  
  if(method == "noise_filter"){
    computedPMDs <- .computePMDsNoiseFilter(methylationData = methylationData, 
                                            regions = regions, 
                                            context = context,
                                            windowSize = windowSize,
                                            kernelFunction = kernelFunction, 
                                            lambda = lambda,
                                            minCytosinesCount = minCytosinesCount,
                                            minMethylation = minMethylation,
                                            maxMethylation = maxMethylation,
                                            minGap = minGap, 
                                            minSize = minSize, 
                                            minReadsPerCytosine = minReadsPerCytosine, 
                                            cores = cores,
                                            BPPARAM = BPPARAM)
  } else if(method == "neighbourhood"){
    computedPMDs <- .computePMDsNeighbourhood(methylationData = methylationData, 
                                              regions = regions, 
                                              context = context, 
                                              minCytosinesCount = minCytosinesCount, 
                                              minMethylation = minMethylation,
                                              maxMethylation = maxMethylation,
                                              minGap = minGap, 
                                              minSize = minSize, 
                                              minReadsPerCytosine = minReadsPerCytosine, 
                                              cores = cores,
                                              BPPARAM = BPPARAM)
  } else if(method == "bins"){
    computedPMDs <- .computePMDsBins(methylationData = methylationData, 
                                     regions = regions, 
                                     context = context, 
                                     binSize = binSize,
                                     minCytosinesCount = minCytosinesCount, 
                                     minMethylation = minMethylation,
                                     maxMethylation = maxMethylation,
                                     minGap = minGap, 
                                     minSize = minSize, 
                                     minReadsPerCytosine = minReadsPerCytosine, 
                                     cores = cores,
                                     BPPARAM = BPPARAM)
    
  } else{
    cat("Unknown method: ",method," \n")
  } 
  
  return(computedPMDs)  
}

#' This function computes the partially methylated regions 
#' using the noise filter method.  
.computePMDsNoiseFilter <- function(methylationData, 
                                    regions = NULL, 
                                    context = "CG",
                                    windowSize = 100,
                                    kernelFunction="triangular", 
                                    lambda=0.5,                                    
                                    minCytosinesCount = 4, 
                                    minMethylation = 0.4, 
                                    maxMethylation = 0.6,
                                    minGap = 200, 
                                    minSize = 50, 
                                    minReadsPerCytosine = 4, 
                                    cores = 1,
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
  .computePMDsInterpolationLoop = function(i){
    computedPMDs <- GRanges()    
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]
      
      cat("Computing PMDs at ",.printGenomicRanges(currentRegion),"\n")
      
      # Select the points in methylationData that we're interested in. These are the 
      # points that lie within 'regions', as well as any that lie within 
      # window.size of them. 
      windowSizeHalf <- floor((windowSize - 1)/2)
      extendedRegion <- currentRegion
      start(extendedRegion) <- start(currentRegion) - windowSizeHalf
      end(extendedRegion) <- end(currentRegion) + windowSizeHalf
      
      overlaps <- findOverlaps(localContextMethylationData, extendedRegion)
      if(length(overlaps) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlaps)]
        
        cat("Calculating interpolations...\n")
        
        #Rcpp
        movingAverageMethylReads <- round(.movingAverage(start(currentRegion), 
                                                          end(currentRegion), 
                                                          start(localMethylationData), 
                                                          localMethylationData$readsM, 
                                                          windowSizeHalf = windowSizeHalf))
        movingAverageTotalReads <- round(.movingAverage(start(currentRegion), 
                                                         end(currentRegion), 
                                                         start(localMethylationData), 
                                                         localMethylationData$readsN, 
                                                         windowSizeHalf = windowSizeHalf))
        movingAverageProportion <- movingAverageMethylReads / movingAverageTotalReads
        
        
        cat("Identifying PMDs...\n")  
        # filtering the expected PMD window
        
        
        # compute the partially methylated cytosines
 
        bufferIndex <- 
          movingAverageTotalReads >= minReadsPerCytosine &
          movingAverageProportion > minMethylation &
          movingAverageProportion < maxMethylation  
        
        bufferIndex[is.na(bufferIndex)] <- FALSE
        PMPs <- sign(bufferIndex)
        
        #join the partially methylated cytosines into regions
        rle <- rle(PMPs)
        rle$cumulative <- cumsum(rle$lengths)
        endOfRuns <- rle$cumulative + start(currentRegion) - 1
        
        PMDs <- GRanges(
          seqnames    = seqnames(currentRegion),
          ranges      = IRanges(endOfRuns - rle$lengths + 1, endOfRuns),
          strand      = strand(currentRegion),
          direction   = rle$values,
          context     = paste(context, collapse = "_")
        )
        
        # Select the crude list of PMDs
        PMDs <- PMDs[!is.na(PMDs$direction) & (PMDs$direction == 1)]
        
        # append current PMDs to the global list of PMDs 
        if(length(computedPMDs) == 0){
          computedPMDs <- PMDs
        } else{
          computedPMDs <- c(computedPMDs,PMDs)
        }
        
      }
    }
    
    return(computedPMDs)
    
  }
  
  # compute the PMDs
  if(cores > 1){
    cat("Compute the PMDs using ", cores, "cores\n")
    computedPMDs <- BiocParallel::bplapply(1:length(regionsList), .computePMDsInterpolationLoop, BPPARAM = BPPARAM)
  } else {
    computedPMDs <- lapply(1:length(regionsList), .computePMDsInterpolationLoop)
  }
  
  
  computedPMDs <- unlist(GRangesList(computedPMDs))
  
  
  if(length(computedPMDs) > 0){
    computedPMDs <- computedPMDs[order(computedPMDs)]
    cat("Analysed reads inside PMDs\n")
    overlaps <- countOverlaps(localContextMethylationData, computedPMDs)
    localContextMethylationDataPMDs <- localContextMethylationData[overlaps > 0]
    if(cores > 1){
      computedPMDsList <- S4Vectors::splitAsList(computedPMDs,  rep(1:cores, length.out=length(computedPMDs)))
      buffercomputedPMDsList <- BiocParallel::bplapply(1:length(computedPMDsList), function(i){ 
        .analyseReadsInsideRegionsPMDs(localContextMethylationDataPMDs, computedPMDsList[[i]])}, 
        BPPARAM = BPPARAM)
      computedPMDs <- unlist(GRangesList(buffercomputedPMDsList))
      computedPMDs <- computedPMDs[order(computedPMDs)]
      
    } else{
      computedPMDs <- .analyseReadsInsideRegionsPMDs(localContextMethylationDataPMDs, computedPMDs)
    }
    computedPMDs <- computedPMDs[!is.na(computedPMDs$cytosinesCount) & computedPMDs$cytosinesCount > 0]

    computedPMDs <- computedPMDs[!is.na(computedPMDs$proportion) &
                                 computedPMDs$proportion >= minMethylation & 
                                   computedPMDs$proportion <= maxMethylation]
    
    cat("Merge PMDs iteratively\n")    
    # Get rid of small gaps between PMDs.
    if(minGap > 0){
      computedPMDs <- .smartMergePMDs(computedPMDs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minReadsPerCytosine = minReadsPerCytosine,
                                      minMethylation = minMethylation,
                                      maxMethylation = maxMethylation,
                                      cores = cores,
                                      BPPARAM = BPPARAM)
    }  
    cat("Filter PMDs \n")    
    if(length(computedPMDs) > 0){
      #remove small PMDs 
      computedPMDs <- computedPMDs[width(computedPMDs) >= minSize]
      if(length(computedPMDs) > 0){
        #remove PMDswith few cytosines
        computedPMDs <- computedPMDs[!is.na(computedPMDs$cytosinesCount) & computedPMDs$cytosinesCount >= minCytosinesCount]
      }
    }
  } 
  if (length(computedPMDs) >0 & length(computedPMDs$direction) != 0){
    computedPMDs$direction <- NULL
  }
  return(computedPMDs)  
}

#' This function computes the partially methylated domains
#' using the neighbourhood method. This assumes the computation of partially methylated 
#' cytosines followed by smart merging of these cytosines while keeping the new PMDs 
#' statistically significant.   
.computePMDsNeighbourhood <- function(methylationData, 
                                      regions = NULL, 
                                      context = "CG",
                                      minCytosinesCount = 4, 
                                      minMethylation = 0.4, 
                                      maxMethylation = 0.6,
                                      minGap = 200, 
                                      minSize = 50, 
                                      minReadsPerCytosine = 4, 
                                      cores = 1,
                                      BPPARAM = BPPARAM) {  
  
  regions <- reduce(regions)
  
  # extract the methylation data in the correct context
  cat("Extract methylation in the corresponding context \n")
  
  
  contextMethylationData <- methylationData[methylationData$context%in%context]
  rm(methylationData)
  localContextMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, regions))]
  rm(contextMethylationData)  
  
  cat("Computing PMDs \n")
  PMPs <- GRanges()
  if(length(localContextMethylationData) > 0){
    PMPs <- localContextMethylationData
    
    PMPs$sumReadsM <- PMPs$readsM
    PMPs$sumReadsN <- PMPs$readsN
    PMPs$proportion <- PMPs$readsM / PMPs$readsN
    PMPs$cytosinesCount <- 1
    PMPs$direction <- sign(PMPs$proportion >= minMethylation &
                             PMPs$proportion <= maxMethylation)
    
    bufferIndex <- 
      PMPs$sumReadsN >=minReadsPerCytosine & PMPs$proportion >= minMethylation &
      PMPs$proportion <= maxMethylation
    PMPs <- PMPs[bufferIndex]
    strand(PMPs) <- "*" 
  }
  
  computedPMDs <- GRanges()
  if(length(PMPs) > 0){    
    cat("Merge PMDs iteratively\n")    
    # Get rid of small gaps between PMDs.
    if(minGap > 0){
      computedPMDs <- .smartMergePMDs(PMPs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minReadsPerCytosine = minReadsPerCytosine, 
                                      minMethylation = minMethylation,
                                      maxMethylation = maxMethylation,
                                      cores = cores,
                                      BPPARAM = BPPARAM)
    } else{
      computedPMDs <- PMPs
    } 
    computedPMDs$ONT_Cm = NULL
    computedPMDs$ONT_C  = NULL
    cat("Filter PMDs \n")    
    if(length(computedPMDs) > 0){
      #remove small PMDs 
      computedPMDs <- computedPMDs[width(computedPMDs) >= minSize]
      if(length(computedPMDs) > 0){
        #remove PMDs with few cytosines
        computedPMDs <- computedPMDs[!is.na(computedPMDs$cytosinesCount) & computedPMDs$cytosinesCount >= minCytosinesCount]
        
      }
      
    }
  }  
  if (length(computedPMDs) >0 & length(computedPMDs$direction) != 0){
    computedPMDs$direction <- NULL
  }
  return(computedPMDs)  
}


#' This function computes the partially methylated regions 
#' using the bins method. 
.computePMDsBins <- function(methylationData,
                             regions = NULL, 
                             context = "CG", 
                             binSize = 100,
                             minCytosinesCount = 4,
                             minMethylation = 0.4, 
                             maxMethylation = 0.6,
                             minGap = 200, 
                             minSize = 50, 
                             minReadsPerCytosine = 4, 
                             cores = 1,
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
  .computePMDsBinsLoop = function(i){
    computedPMDs <- GRanges()
    for(index in 1:length(regionsList[[i]])){
      currentRegion <- regionsList[[i]][index]
      
      
      cat("Computing PMDs at ",.printGenomicRanges(currentRegion),"\n")
      
      seqs <- seq(start(currentRegion), (end(currentRegion)-binSize), by = binSize);
      
      bins <- GRanges(seqnames(currentRegion), IRanges(seqs, (seqs+binSize-1)))
      
      overlapsBins <- findOverlaps(localContextMethylationData, currentRegion)
      
      if(length(overlapsBins) > 0){
        localMethylationData <- localContextMethylationData[queryHits(overlapsBins)]
        
        cat("Count inside each bin...\n")
        #bins <- .analyseReadsInsideRegions(localMethylationData, bins, context, cores)
        bins <- .analyseReadsInsideBinsPMDs(localMethylationData, bins, currentRegion)
        
        cat("Filter the bins...\n")
        # Get rid of the bins with fewer than minCytosinesCount cytosines inside.  
        bins  <- bins[bins$cytosinesCount >= minCytosinesCount]
        
        # Get rid of the bins with fewer than minReadsPerCytosine reads per cytosine.  
        bins  <- bins[(bins$sumReadsN/bins$cytosinesCount >= minReadsPerCytosine)]
        
        # Get rid of the bins with outlier about ranges for PMD parameter
        bins  <- bins[bins$proportion >= minMethylation & bins$proportion <= maxMethylation]
        
        
        cat("Identifying PMDs...\n")    
        bins$context <- rep(paste(context, collapse = "_"), length(bins))
        bins$direction <- rep(NA, length(bins))
        bins$direction <- sign(bins$proportion)
        
        # Select the crude list of PMDs
        PMDs <- bins[!is.na(bins$direction) & bins$direction == 1]
        
        # append current PMDs to the global list of PMDs 
        if(length(computedPMDs) == 0){
          computedPMDs <- PMDs
        } else{
          computedPMDs <- c(computedPMDs,PMDs)
        }
      }
    }
    return(computedPMDs)
  }
  
  # compute the PMDs
  if(cores > 1){
    cat("Compute the PMDs using ", cores, "cores\n")
    computedPMDs <- BiocParallel::bplapply(1:length(regionsList), .computePMDsBinsLoop, BPPARAM = BPPARAM)
  } else {
    computedPMDs <- lapply(1:length(regionsList), .computePMDsBinsLoop)
  }
  
  
  computedPMDs <- unlist(GRangesList(computedPMDs))
  
  if(length(computedPMDs) > 0){
    
    cat("Merge adjacent PMDs\n")    
    computedPMDs <- computedPMDs[order(computedPMDs)]
   
    cat("Merge PMDs iteratively\n")    
    # Get rid of small gaps between PMDs.
    if(minGap > 0){
      computedPMDs <- .smartMergePMDs(computedPMDs, 
                                      minGap = minGap, 
                                      respectSigns = TRUE, 
                                      methylationData = localContextMethylationData,
                                      minReadsPerCytosine = minReadsPerCytosine,
                                      minMethylation = minMethylation,
                                      maxMethylation = maxMethylation,
                                      cores = cores,
                                      BPPARAM = BPPARAM)
    }  
    
    cat("Filter PMDs \n")    
    if(length(computedPMDs) > 0){
      #remove small PMDs 
      computedPMDs <- computedPMDs[width(computedPMDs) >= minSize]
      if(length(computedPMDs) > 0){
        #remove PMDswith few cytosines
        computedPMDs <- computedPMDs[!is.na(computedPMDs$cytosinesCount) & computedPMDs$cytosinesCount >= minCytosinesCount]
      }
    }
  }  
  if (length(computedPMDs) >0 & length(computedPMDs$direction) != 0){
    computedPMDs$direction <- NULL
  }
  return(computedPMDs)
}

#' This function verifies whether a set of potential PMDs (e.g. genes, 
#' transposons, CpG islands) are partially methylated or not.
#'
#' @title Filter PMDs 
#' @param methylationData the methylation data in condition
#' (see \code{\link{ontSampleGRangesList}}).
#' @param potentialPMDs a \code{\link[GenomicRanges]{GRanges}} object with potential PMDs 
#' where to compute the PMDs. This can be a a list of gene and/or transposable 
#' elements coordinates.
#' @param context the context in which the PMDs are computed (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param minCytosinesCount PMDs with less cytosines in the specified context 
#' than \code{minCytosinesCount} will be discarded.
#' @param minMethylation Numeric [0,1]; minimum mean methylation proportion.
#' @param maxMethylation Numeric [0,1]; maximum mean methylation proportion.
#' @param minReadsPerCytosine  PMDs with the average number of reads lower than 
#' \code{minReadsPerCytosine} are discarded. 
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return a \code{\link[GenomicRanges]{GRanges}} object with 5 metadata columns that contain 
#' the PMDs; see \code{\link{computePMDs}}.
#' @seealso \code{\link{PMDsNoiseFilterCG}}, \code{\link{computePMDs}}, 
#' \code{\link{analyseReadsInsideRegionsForCondition}}  
#' and \code{\link{mergePMDsIteratively}}
#' @examples
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
#' # filter genes that are partially methylated in the two conditions
#' PMDsGenesCG <- filterPMDs(ontSampleGRangesList[["GM18501"]], 
#'                potentialPMDs = transcript, 
#'                context = "CG", minMethylation = 0.4, maxMethylation = 0.6,
#'                minCytosinesCount = 4, minReadsPerCytosine = 3, cores = 1)
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @export
filterPMDs <- function(methylationData, 
                       potentialPMDs, 
                       context = "CG",   
                       minCytosinesCount = 4, 
                       minMethylation = 0.4, 
                       maxMethylation = 0.6,
                       minReadsPerCytosine = 3, 
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
  
  .validateMethylationData(methylationData, variableName="methylationData")
  
  regions <- getWholeChromosomes(methylationData)
  
  .validateContext(context)
  
  .validateGRanges(potentialPMDs, generateGenomeWide=FALSE, variableName="potentialPMDs", minLength=NULL)
  
  regions <- .validateGRanges(regions, methylationData)
  
  .stopIfNotAll(c(.isInteger(minCytosinesCount, positive=TRUE)), 
                " the minCytosinesCount is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum number of reads in a bin is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(minMethylation), is.numeric(minMethylation), minMethylation >= 0), 
                " minMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(maxMethylation), is.numeric(maxMethylation), maxMethylation >= 0), 
                " maxMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(maxMethylation >= minMethylation), 
                " maxMethylation should be higher than minMethylation value")
  
  regions <- reduce(regions)
  
  if(length(potentialPMDs) > 0){
    
    
    
    # extract the methylation data in the correct context
    cat("Extract methylation in the corresponding context \n")
    
    contextMethylationData <- methylationData[methylationData$context%in%context]
    rm(methylationData)
    localContextMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, regions))]
    rm(contextMethylationData)  
    
    regionsList <- .splitGRangesEqualy(regions, cores)
    
    # inner loop function for BiocParallel::bplapply
    .filterPMDsLoop = function(i){
      computedPMDs <- GRanges()  
      for(index in 1:length(regionsList[[i]])){
        currentRegion <- regionsList[[i]][index]
        
        
        cat("Computing PMDs at ",.printGenomicRanges(currentRegion),"\n")
        
        cat("Selecting data...\n")
        
        # Select the points in methylationData that we're interested in. These are the 
        # points that lie within 'regions', as well as any that lie within 
        # window.size of them. 
        
        overlapsPotentialPMDs <- findOverlaps(potentialPMDs, currentRegion)
        if(length(overlapsPotentialPMDs) > 0){
          potentialPMDsLocal <- potentialPMDs[queryHits(overlapsPotentialPMDs)]
          
          localMethylationData <- localContextMethylationData[queryHits(findOverlaps(localContextMethylationData, currentRegion))]
          potentialPMDsLocal <- .analyseReadsInsideRegionsPMDs(localMethylationData, potentialPMDsLocal)
          
          if(length(computedPMDs) == 0){
            computedPMDs <- potentialPMDsLocal
          } else{
            computedPMDs <- c(computedPMDs,potentialPMDsLocal)
          }
        } 
      }
      return(computedPMDs)
    }
    
    # compute the PMDs
    if(cores > 1){
      cat("Compute the PMDs using ", cores, "cores\n")
      computedPMDs <- BiocParallel::bplapply(1:length(regionsList), .filterPMDsLoop, BPPARAM = BPPARAM)
    } else {
      computedPMDs <- lapply(1:length(regionsList), .filterPMDsLoop)
    }
    
    
    computedPMDs <-  unlist(GRangesList(computedPMDs))
    
    if(length(computedPMDs) > 0){
      cat("Identifying PMDs...\n")    
      
      bufferIndex <- !is.na(computedPMDs$proportion) &
        computedPMDs$proportion >= minMethylation & 
        computedPMDs$proportion <= maxMethylation &
        computedPMDs$sumReadsN/computedPMDs$cytosinesCount >= minReadsPerCytosine &
        computedPMDs$cytosinesCount >= minCytosinesCount
      
      computedPMDs <- computedPMDs[bufferIndex]  
    }
  } else{
    computedPMDs <- GRanges() 
  }
  
  if(length(computedPMDs) > 0){
    computedPMDs <- computedPMDs[order(computedPMDs)]
  }
  
  return(computedPMDs)
  
}



#' This function takes a list of PMDs and attempts to merge PMDs while keeping 
#' the new PMDs statistically significant.
#'
#' @title Merge PMDs iteratively
#' @param PMDs the list of PMDs as a \code{\link[GenomicRanges]{GRanges}} object; 
#' e.g. see \code{\link{computePMDs}}
#' @param minGap PMDs separated by a gap of at least \code{minGap} are not 
#' merged.
#' @param respectSigns logical value indicating whether to respect the sign when 
#' joining PMDs.
#' @param methylationData the methylation data in GRanges
#' (see \code{\link{ontSampleGRangesList}}).
#' @param context the context in which the PMDs are computed (\code{"CG"}, 
#' \code{"CHG"} 
#' or \code{"CHH"}).
#' @param minReadsPerCytosine two adjacent PMDs are merged only if the number of 
#' reads per cytosine of the new DMR is higher than \code{minReadsPerCytosine}.
#' @param minMethylation Numeric [0,1]; minimum mean methylation proportion.
#' @param maxMethylation Numeric [0,1]; maximum mean methylation proportion.
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return the reduced list of PMDs as a \code{\link[GenomicRanges]{GRanges}} object; 
#' e.g. see \code{\link{computePMDs}}
#' @seealso \code{\link{filterPMDs}}, \code{\link{computePMDs}}, 
#' \code{\link{analyseReadsInsideRegionsForCondition}} and 
#' \code{\link{PMDsNoiseFilterCG}}
#' @examples
#' # load the ONT methylation data 
#' data(ontSampleGRangesList)
#' 
#' # load the PMDs in CG context they were computed with minGap = 200
#' data(PMDsNoiseFilterCG)
#' 
#'
#' # merge the PMDs 
#' PMDsNoiseFilterCGLarger <- mergePMDsIteratively(PMDsNoiseFilterCG[1:100], 
#'                            minGap = 500, respectSigns = TRUE, 
#'                            ontSampleGRangesList[["GM18501"]], context = "CG", 
#'                            minReadsPerCytosine = 4, minMethylation = 0.4, 
#'                            maxMethylation = 0.6, cores = 1)
#' 
#' \dontrun{
#' # set genomic coordinates where to compute PMDs
#' chr1_ranges <- GRanges(seqnames = Rle("chr1"), ranges = IRanges(1E6+5E5,2E6))
#' 
#' # compute PMDs and remove gaps smaller than 200 bp
#' PMDsNoiseFilterCG200 <- computePMDs(ontSampleGRangesList[["GM18501"]], 
#'                        regions = chr1_ranges, context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular",  
#'                        minCytosinesCount = 1, minMethylation = 0.4, 
#'                        maxMethylation = 0.6, minGap = 0, minSize = 200, 
#'                        minReadsPerCytosine = 1, cores = 1)
#' PMDsNoiseFilterCG0 <- computePMDs(ontSampleGRangesList[["GM18501"]], 
#'                        regions = chr1_ranges, context = "CG", method = "noise_filter", 
#'                        windowSize = 100, kernelFunction = "triangular", 
#'                        minCytosinesCount = 1, minMethylation = 0.4, 
#'                        maxMethylation = 0.6, minGap = 0, minSize = 0, 
#'                        minReadsPerCytosine = 1, cores = 1)
#' PMDsNoiseFilterCG0Merged200 <- mergePMDsIteratively(PMDsNoiseFilterCG0, 
#'                              minGap = 200, respectSigns = TRUE, 
#'                              ontSampleGRangesList[["GM18501"]], context = "CG",  
#'                              minReadsPerCytosine = 4, minMethylation = 0.4, 
#'                              maxMethylation = 0.6, cores = 1)                      
#' 
#' #check that all newley computed PMDs are identical
#' print(all(PMDsNoiseFilterCG200 == PMDsNoiseFilterCG0Merged200))
#' 
#' }
#' 
#' @author Nicolae Radu Zabet and Young Jun Kim
#' 
#' @export
mergePMDsIteratively <- function(PMDs, 
                                 minGap = 200, 
                                 respectSigns = TRUE, 
                                 methylationData,
                                 context = "CG",
                                 minReadsPerCytosine = 4, 
                                 minMethylation = 0.4, 
                                 maxMethylation = 0.6,
                                 parallel = FALSE,
                                 BPPARAM = NULL,
                                 cores = NULL){
  
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
  
  .validateMethylationData(methylationData, variableName="methylationData")
  .validateContext(context)
  .validateGRanges(PMDs, generateGenomeWide=FALSE, variableName="PMDs", minLength=NULL)

  .stopIfNotAll(c(.isInteger(minGap, positive=TRUE)),
                " the minimum gap between PMDs is an integer higher or equal to 0")
  
  .stopIfNotAll(c(.isInteger(minReadsPerCytosine, positive=TRUE)), 
                " the minimum average number of reads in a DMR is an integer higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(minMethylation), is.numeric(minMethylation), minMethylation >= 0), 
                " minMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(!is.null(maxMethylation), is.numeric(maxMethylation), maxMethylation >= 0), 
                " maxMethylation needs to be a numeric value higher or equal to 0")
  
  .stopIfNotAll(c(maxMethylation >= minMethylation), 
                " maxMethylation should be higher than minMethylation value")
  
  totalRegion <- reduce(PMDs, drop.empty.ranges=FALSE, min.gapwidth=minGap, ignore.strand=TRUE)
  
  contextMethylationData <- methylationData[methylationData$context%in%context]
  rm(methylationData)
  localContextMethylationData <- contextMethylationData[queryHits(findOverlaps(contextMethylationData, totalRegion))]
  rm(contextMethylationData)  
  
  cat("Merge PMDs iteratively ...\n")
  
  
  return(.smartMergePMDs(PMDs, 
                         minGap = minGap, 
                         respectSigns = respectSigns, 
                         methylationData = localContextMethylationData,
                         minReadsPerCytosine = minReadsPerCytosine,
                         minMethylation = minMethylation,
                         maxMethylation = maxMethylation,
                         cores = cores,
                         BPPARAM = BPPARAM))                                               
}


#' This function extracts from the methylation data the total number of reads, 
#' the number of methylated reads and the number of cytosines in the specific 
#' context from a region (e.g. PMDs)
#'
#' @title Analyse reads inside regions for condition
#' @param regions a \code{\link[GenomicRanges]{GRanges}} object with a list of regions on the 
#' genome; e.g. could be a list of PMDs
#' @param methylationData the methylation data in one condition
#' (see \code{\link{ontSampleGRangesList}}).
#' @param context the context in which to extract the reads (\code{"CG"}, 
#' \code{"CHG"} or \code{"CHH"}).
#' @param label a string to be added to the columns to identify the condition
#' @param parallel Logical; run in parallel if \code{TRUE}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallel execution.
#'    This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#' @param cores Integer number of workers (must not exceed BPPARAM$workers).
#'    This value will automatically set as the maximum number of system workers,
#'    also able to set as manually.
#' @return a \code{\link[GenomicRanges]{GRanges}} object with additional four metadata columns
#' \describe{
#'  \item{sumReadsM}{the number of methylated reads}
#'  \item{sumReadsN}{the total number of reads} 
#'  \item{proportion}{the proportion methylated reads} 
#'  \item{cytosinesCount}{the number of cytosines in the regions} 
#' }
#' @seealso \code{\link{filterPMDs}}, \code{\link{computePMDs}}, 
#' \code{\link{PMDsNoiseFilterCG}}, and \code{\link{mergePMDsIteratively}}
#' @examples
#' 
#' # load the ONT methylation data
#' data(ontSampleGRangesList)
#'  
#' #load the PMDs in CG context. These PMDs were computed with minGap = 200.
#' data(PMDsNoiseFilterCG)
#' 
#' #retrive the number of reads in CG context in GM18501
#' PMDsNoiseFilterCGreadsCG <- analyseReadsInsideRegionsForConditionPMD(
#'                              PMDsNoiseFilterCG[1:10], 
#'                              ontSampleGRangesList[["GM18501"]], context = "CG", 
#'                              label = "GM18501")
#' 
#' 
#' @author Nicolae Radu Zabet and Young Jun Kim
#' 
#' @export
analyseReadsInsideRegionsForConditionPMD <- function(regions,
                                                  methylationData,
                                                  context,
                                                  label = "",
                                                  parallel = FALSE,
                                                  BPPARAM = NULL,
                                                  cores = NULL){
  
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
  
  .validateGRanges(regions, generateGenomeWide=FALSE, variableName="regions", minLength=NULL)
  .validateMethylationData(methylationData, variableName="methylationData")
  .validateContext(context)
  
  cat("Extract methylation levels in corresponding context ...\n")
  contextMethylationData <- methylationData[methylationData$context%in%context]
  rm(methylationData)
  
  cat("Compute reads inside each region ...\n")
  if(length(regions) > 0){
    if(cores > 1){
      cat("Analyse reads in regions using ", cores, "cores\n")
      regionsList <- split(regions, rep(1:cores, length.out=length(regions)))
      .analyseReadsInsideRegionsForConditionLoop = function(i){
        regionsLocal <- .analyseReadsInsideRegionsForCondition(regionsList[[i]], 
                                                               contextMethylationData, 
                                                               label = label,
                                                               context = context)
        return(regionsLocal)
      }
      regions <- BiocParallel::bplapply(1:length(regionsList), 
                                    .analyseReadsInsideRegionsForConditionLoop, BPPARAM = BPPARAM)
      regions <- unlist(GRangesList(regions))
    } else{
      regions <- .analyseReadsInsideRegionsForCondition(regions, 
                                                        contextMethylationData, 
                                                        label = label,
                                                        context = context)
    }
    regions <- regions[order(regions)]
  }
  return(regions)
}  

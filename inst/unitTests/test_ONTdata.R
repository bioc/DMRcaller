library(DMRcaller)
library(RUnit)

test_readONTbam_chr1_sample <- function() {
  ## 1) Locate the example BAM and its index in extdata/
  bamfile <- system.file("extdata", "scanBamChr1Random5.bam",
                         package = "DMRcaller")
  stopifnot(file.exists(bamfile),
            file.exists(paste0(bamfile, ".bai")))
  
  
  ## 2) Call readONTbam with ref_gr = NULL so it builds its own
  ##    cytosine reference for CG sites on chr1
  gr <- readONTbam(
    bamfile = bamfile,
    ref_gr  = NULL,
    context = "CG",
    chr     = "chr1"
  )
  
  ## test1. Sanity check: total coverage should be at least as large
  ##    as the count of modified reads at every site
  checkTrue(
    all(gr$readsN >= gr$readsM)
  )
  
  ## test2. selectCytosine function check (reference GRanges)
  checkTrue(unique(seqnames(gr))=="chr1")
  checkTrue(unique(gr$context)=="CG")
}


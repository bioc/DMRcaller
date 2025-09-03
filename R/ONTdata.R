#' This function enumerate cytosine positions in a BSgenome reference
#'
#' @title Select Cytosine Positions
#' @description
#' Constructs a \code{\link[GenomicRanges]{GRanges}} of all cytosine positions
#' in the specified BSgenome (or BSgenome package name), optionally filtering
#' by methylation context (\code{"CG"}, \code{"CHG"}, \code{"CHH"}), by
#' chromosome, and by genomic region.
#'
#' @param genome
#'   A \code{BSgenome} object or the name (character) of a BSgenome package
#'   to use as the reference genome. If a package name is given, it will be loaded
#'   automatically (default: \code{BSgenome.Hsapiens.UCSC.hg38}). 
#'   \strong{Note:} When running on an HPC system, please ensure that the required
#'   \code{BSgenome} package is already installed and loaded in advance.
#' @param context
#'   A character vector of one or more methylation contexts to include:
#'   \code{"CG"}, \code{"CHG"}, and/or \code{"CHH"}.  Defaults to all three.
#' @param chr
#'   An optional character vector of chromosome names to restrict the
#'   enumeration.  If \code{NULL}, all sequences in \code{genome} are used.
#' @param region
#'   An optional \code{\link[GenomicRanges]{GRanges}} object specifying
#'   subregions to keep.  Requires \code{chr} to be non-\code{NULL}.
#'
#' @return
#' A \code{\link[GenomicRanges]{GRanges}} object with one 1-bp range per
#' cytosine, and two metadata columns:
#' \describe{
#'   \item{context}{A factor indicating the context (\code{"CG"},
#'     \code{"CHG"}, or \code{"CHH"}).}
#'   \item{trinucleotide_context}{Character or factor giving the
#'     surrounding trinucleotide sequence (or \code{NA} for \code{"CG"}).}
#' }
#' @seealso \code{\link{readONTbam}}, \code{\link{ontSampleGRangesList}} 
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # All cytosines in hg38:
#' gr_all <- selectCytosine()
#'
#' # Only CpG sites on chr1 and chr2:
#' gr_chr1_2 <- selectCytosine(context="CG", chr=c("chr1","chr2"))
#'
#' # CHH sites in a specific region on chr3:
#' my_region <- GRanges("chr3", IRanges(1e6, 1e6 + 1e5))
#' gr_region <- selectCytosine(context="CHH", chr="chr3", region=my_region)
#' }
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @import GenomicRanges Biostrings
#' @export
selectCytosine <- function(genome   = BSgenome.Hsapiens.UCSC.hg38,
                           context  = c("CG","CHG","CHH"),
                           chr      = NULL,
                           region   = NULL) {
  ##Parameters checking
  cat("Parameters checking ...\n")
  
  .validateContext(context)
  
  .validateGenome(genome)
  
  .validateChromosome(chr, genome)
  
  ## select the cytosine based on the reference genome
  cat("Select the Cytosine Positions in the specified BSgenome ...\n")
  ## --- literal 2-/3-mer definitions
  pat_dict <- list(
    CG  = "CG",
    CHG = c("CAG","CCG","CTG"),
    CHH = c("CAA","CAC","CAT",
            "CCA","CCC","CCT",
            "CTA","CTC","CTT")
  )
  if (!all(context %in% names(pat_dict))) {
    stop("Unknown context: ", paste(setdiff(context, names(pat_dict)), collapse=", "))
  }
  
  chroms <- if (is.null(chr)) seqnames(genome) else intersect(chr, seqnames(genome))
  if (!length(chroms)) stop("No valid chromosomes selected.")
  
  per_chr <- lapply(chroms, function(chrN) {
    seqN <- genome[[chrN]]
    seql <- seqlengths(genome)[chrN]
    ctx_hits <- lapply(intersect(context, names(pat_dict)), function(ctx) {
      patterns <- pat_dict[[ctx]]
      
      # Forward strand matches
      pos_fwd <- unique(unlist(lapply(patterns, function(pat) start(matchPattern(pat, seqN))), use.names=FALSE))
      pos_fwd <- pos_fwd[pos_fwd >= 1 & pos_fwd <= (seql - ifelse(ctx=="CG",1,2))]
      gr_fwd <- if (length(pos_fwd)) {
        gr <- GRanges(chrN, IRanges(start=pos_fwd, width=1L), strand="+")
        tri_fwd <- if (ctx=="CG") rep(NA_character_, length(pos_fwd))
        else as.character(getSeq(genome, chrN, pos_fwd, width=3L))
        mcols(gr) <- DataFrame(context=factor(ctx, levels=names(pat_dict)),
                               trinucleotide_context=tri_fwd)
        gr
      } else GRanges()
      
      # Reverse strand for CG/CHG
      gr_rev <- GRanges()
      if (ctx %in% c("CG","CHG") && length(gr_fwd)) {
        shift <- if (ctx=="CG") 1L else 2L
        gr <- GRanges(seqnames(gr_fwd), IRanges(start(gr_fwd)+shift, width=1L), strand="-")
        tri_rev <- if (ctx=="CG") rep(NA_character_, length(gr_fwd))
        else as.character(reverseComplement(DNAStringSet(mcols(gr_fwd)$trinucleotide_context)))
        mcols(gr) <- DataFrame(context=mcols(gr_fwd)$context,
                               trinucleotide_context=tri_rev)
        if (ctx %in% "CHG"){
          valid <- mcols(gr)$trinucleotide_context %in% pat_dict[["CHG"]]
          gr<- gr[valid]
        }
        gr_rev <- gr
      }
      
      # Reverse strand for CHH
      gr_rev_chh <- GRanges()
      if (ctx=="CHH") {
        patterns_rev <- as.character(reverseComplement(DNAStringSet(patterns)))
        pos_rev <- unique(unlist(lapply(patterns_rev, function(pat) start(matchPattern(pat, seqN))), use.names=FALSE))
        pos_rev <- pos_rev[pos_rev >=1 & pos_rev <= (seql-2)]
        if (length(pos_rev)) {
          gr2 <- GRanges(chrN, IRanges(start=pos_rev+2, width=1L), strand="-")
          tri_fwd <- as.character(getSeq(genome, chrN, pos_rev, width=3L))
          tri_rev <- as.character(reverseComplement(DNAStringSet(tri_fwd)))
          mcols(gr2) <- DataFrame(context=factor(ctx, levels=names(pat_dict)),
                                  trinucleotide_context=tri_rev)
          gr_rev_chh <- gr2
        }
      }
      
      c(gr_fwd, gr_rev, gr_rev_chh)
    })
    do.call(c, ctx_hits)
  })
  
  # Harmonize seqlevels
  per_chr <- lapply(per_chr, function(gr) {
    seqlevels(gr)  <- chroms
    seqlengths(gr) <- seqlengths(genome)[chroms]
    gr
  })
  
  # Flatten and sort
  ref_gr <- sort(do.call(c, per_chr))
  o <- order(as.character(seqnames(ref_gr)), start(ref_gr))
  ref_gr <- ref_gr[o]
  
  # If region specified, subset final GRanges
  if (!is.null(region)) {
    # Ensure seqlevels compatibility
    seqlevels(region, pruning.mode="coarse") <- seqlevels(ref_gr)
    ref_gr <- subsetByOverlaps(ref_gr, region, ignore.strand=FALSE)
  }
  
  ref_gr
}

#' This function read and annotate ONT MM/ML tags against a cytosine reference
#'
#' @title  Load ONT BAM, decode MM/ML, and count modified vs. unmodified reads
#' @description
#' \code{readONTbam()} takes an indexed Nanopore BAM file with MM/ML tags,
#' decodes each read’s per-C modification probabilities, and overlays
#' them on a \code{GRanges} of candidate cytosine sites.
#' It returns a copy of \code{ref_gr} augmented with:
#' \itemize{
#'   \item \code{ONT_Cm} — comma-delimited read‐indices called “modified”  
#'   \item \code{ONT_C}  — comma-delimited read‐indices covering but _not_ modified  
#'   \item \code{readsM} — count of modified reads at each site  
#'   \item \code{readsN} — total same-strand coverage at each site  
#' }
#'   
#' You can either supply your own \code{ref_gr} (e.g.\ from \code{selectCytosine()})
#' or leave it \code{NULL} and pass \code{context}, \code{chr}, \code{region}
#' to build \code{ref_gr} on the fly.
#'
#' @param bamfile      Path to an indexed ONT BAM with MM/ML tags.
#' @param ref_gr       A \code{GRanges} of genomic cytosine positions to annotate.  
#'                     If \code{NULL}, will be created via \code{selectCytosine()}  
#'                     using \code{context}, \code{chr}, \code{region}.  
#' @param modif        Character vector of MM codes to treat as “modified”  
#'                     (e.g. \code{"C+m?"}, \code{"C+h?"}, \code{"C+m."}).  
#' @param prob_thresh  Numeric in \[0,1\] — minimum ML probability to call a read  
#'                     “modified.”  
#' @param genome       A \code{BSgenome} object such as \code{BSgenome.Hsapiens.UCSC.hg38}.
#'                     This is used to extract sequence context and must be loaded in advance.
#'                     \strong{Note:} When running on an HPC system, 
#'                     please ensure that the required
#'                     \code{BSgenome} package is already installed and loaded in advance.  
#' @param context      Sequence context for \code{selectCytosine()}  
#'                     (e.g. \code{"CG"}, \code{"CHG"}, \code{"CHH"}).  
#' @param chr          Chromosome names to restrict \code{selectCytosine()}.  
#' @param region       A \code{GRanges} to further subset \code{selectCytosine()}.
#' @param synonimous   Logical (default: FALSE). If TRUE, include modified calls that match the specified context sequence (e.g. CGG),
#'                     even if the site was previously excluded due to deletion or mismatch.
#'                     For example, if a deletion occurs at position 234523, but the surrounding context still forms CGG,
#'                     then the modified C at 234523 will be retained (Nmod=1).
#' @param parallel     Logical. If TRUE, automatically detect your system condition and
#'                     decoding will use parallel threads via BiocParallel::. 
#'                     If FALSE (default), decoding is done serially.
#' @param BPPARAM      A \code{BiocParallelParam} object controlling parallel execution.
#'                     This value will automatically set when parallel is \code{TRUE}, also able to set as manually.
#'
#'
#' @return A \code{GRanges} of the same length as \code{ref_gr}, with four
#'   additional metadata columns:
#'   \describe{
#'     \item{ONT_Cm}{comma-delimited read‐indices called modified}
#'     \item{ONT_C}{comma-delimited read‐indices covering but unmodified}
#'     \item{readsM}{integer count of modified reads per site}
#'     \item{readsN}{integer count of same‐strand reads covering each site}
#'   }
#' @seealso \code{\link{selectCytosine}}, \code{\link{computeDMRs}}, 
#' \code{\link{computePMDs}}, \code{\link{computeCoMethylation}}, 
#' \code{\link{filterVMRsONT}}, \code{\link{ontSampleGRangesList}},
#' \code{\link{scanBamChr1Random5}} 
#' @examples
#' \dontrun{
#' library(DMRcaller)
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' 
#' # set the bam file directory
#' bam_path <- system.file("extdata", "scanBamChr1Random5.bam", package="DMRcaller")
#' 
#' # read ONTbam file (chromosome 1 only) in CG context with BSgenome.Hsapiens.UCSC.hg38
#' ONTSampleGRanges <- readONTbam(bamfile = bam_path, ref_gr = NULL, modif = "C+m?",
#'                          prob_thresh = 0.50,genome = BSgenome.Hsapiens.UCSC.hg38,
#'                          context = "CG", chr = "chr1", region = NULL,
#'                          synonymous = FALSE, parallel = FALSE, BPPARAM = NULL)
#' }
#'
#' @author Nicolae Radu Zabet and Young Jun Kim
#' @import GenomicAlignments
#' @import Rsamtools
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps mcols
#' @importFrom IRanges IRanges IntegerList
#' @importFrom Biostrings reverseComplement matchPattern getSeq
#' @importFrom BiocParallel SerialParam bplapply
#' @export
readONTbam <- function(bamfile,
                       ref_gr      = NULL,
                       modif       = "C+m?",
                       prob_thresh = 0.50,
                       genome      = BSgenome.Hsapiens.UCSC.hg38,
                       context     = "CG",
                       chr         = NULL,
                       region      = NULL,
                       synonymous  = FALSE, 
                       parallel    = FALSE,
                       BPPARAM     = NULL) {
  ### PUT THE include_diff = FALSE, include_nocall = FALSE) after including parameter 
 
  ##Parameters checking
  cat("Parameters checking ...\n")
  
  .validateContext(context)
  
  .validateBamfile(bamfile)
  
  .validateGenome(genome)
  
  .validateChromosome(chr, genome)

  .validateModif(modif, bamfile, chr, genome)

  if (is.null(ref_gr)){
    ref_gr = selectCytosine(genome, context, chr, region)
  }
  region <- .validateGRanges(region, ref_gr)
 
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
  
  .stopIfNotAll(c(!is.null(prob_thresh), 
                  is.numeric(prob_thresh), 
                  prob_thresh > 0, 
                  prob_thresh < 1),
                " the p-value threshold needs to be in the interval (0,1)")
  
  .stopIfNotAll(c(!is.null(synonymous), 
                  is.logical(synonymous)),
                " the synonymous needs to be logical")
  if (synonymous == TRUE){
    .stopIfNotAll(c(!is.null(synonymous), 
                    is.logical(synonymous),
                    context %in% c("CHG", "CHH")),
                  " the context should choose the CHG or CHH context when synonymous is TRUE \n set the synonymous value as FALSE when context is CG")
  }

  
   # @@@@time t1
  t1 <- proc.time()

  ## 1) preprocess the scanned bam list with filtering 
  # existing indexing + ScanBamParam code
  cat("[readONTbam] Generating index file for BAM ...\n")
  if (!file.exists(paste0(bamfile,".bai")))
    Rsamtools::indexBam(bamfile)
  cat("[readONTbam] Loading BAM ...\n")
  if (!is.null(chr) && length(chr) > 0) {
    which_region <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, seqlengths(genome)[chr]))
    param <- Rsamtools::ScanBamParam(
      what = c("qname", "rname", "pos", "cigar", "flag", "seq", "mapq"),
      tag  = c("MM", "ML"),
      which = which_region
    )
  } else {
    param <- Rsamtools::ScanBamParam(
      what = c("qname", "rname", "pos", "cigar", "flag", "seq", "mapq"),
      tag  = c("MM", "ML")
    )
  }
  sb <- Rsamtools::scanBam(bamfile, param = param)[[1]]
  # @@@@@time t2
  t2 <- proc.time()
  cat(sprintf("Step 1 (scanBam): %.2f sec\n", (t2-t1)[["elapsed"]]))
  
  # Filter out 
  cat("[readONTbam] Filtering invalid or empty reads ...\n")
  valid_cigar <- !is.na(sb$cigar)
  valid_seq   <- nchar(sb$seq) > 0
  # handle primary/secondary/supplementary filtering
  is_secondary      <- bitwAnd(sb$flag, 0x100) != 0
  is_supplementary  <- bitwAnd(sb$flag, 0x800) != 0
  has_tag           <- !vapply(sb$tag$ML, function(x) all(is.na(x)), logical(1))
  has_hardclip <- grepl("H", sb$cigar)
  
  # primary: not secondary and not supplementary
  is_primary <- !is_secondary & !is_supplementary
  
  # keep secondary/supplementary 
  valid_secondary_supplementary <- is_primary | ((is_secondary | is_supplementary) & has_tag & !has_hardclip & valid_seq)
  
  # keep only those passing both filters
  if (!is.null(chr)){
    keep <- valid_cigar & sb$rname %in% chr & valid_secondary_supplementary
  } else {
    keep <- valid_cigar & sb$rname & !is.na(sb$rname) & valid_secondary_supplementary
  }
  keep[is.na(keep)] <- FALSE
  
  # subset all vectors in sb
  for (nm in c("qname","rname","pos","cigar","flag","seq","mapq")) {
    sb[[nm]] <- sb[[nm]][keep]
  }
  # subset all BAM tags
  if (!is.null(sb$tag)) {
    sb$tag <- lapply(sb$tag, `[`, keep)
  }
  
  # Decode the strand from the FLAG bit. 
  # In SAM/BAM, bit 0x10 (decimal 16) is set if the read maps to the reverse strand.
  strands <- ifelse(bitwAnd(sb$flag, 0x10) != 0, "-", "+")
  # 
  # # reconstruct aligned blocks per read
  # cat("[readONTbam] Building coverage blocks per read ...\n")
  # cigar_str <- sb$cigar
  # blocks_ir <- cigarRangesAlongReferenceSpace(
  #   cigar_str,
  #   pos = sb$pos,
  #   ops = c("M", "=", "X")
  # )
  # blocks_grl <- GRangesList(
  #   mapply(function(rn, ir, st) {
  #     GRanges(seqnames = rn,
  #             ranges   = ir,
  #             strand   = st)
  #   },
  #   rn = sb$rname,
  #   ir = blocks_ir,
  #   st = strands,
  #   SIMPLIFY = FALSE)
  # )
  
  # After sb is loaded and filtered
  aln <- GAlignments(seqnames=sb$rname, pos=sb$pos, cigar=sb$cigar, seq=sb$seq)
  ref_seq_list <- getSeq(genome, seqnames(aln), start(aln), end(aln))
  aligned_seq_list <- sequenceLayer(sb$seq, cigar=sb$cigar, from="query", to="reference")
  
  # build a lookup: read name -> its index ; easy to navigating less memory
  read_order <- seq_along(sb$qname)
  names(read_order) <- sb$qname
  # @@@@ time t3
  t3 <- proc.time()
  cat(sprintf("Step 2 (filter): %.2f sec\n", (t3-t2)[["elapsed"]]))
  
  # --------------
  #  Helper: decode one read’s MM/ML using C‐counting
  # --------------
  # *1) Helper: extract the number of soft‐clipped bases from a CIGAR string
  # Vectorized function to extract soft‐ and hard‐clip counts from a CIGAR string
  
  .getClipping <- function(cigar) {
    # Leading clips: soft (S) plus hard (H)
    lead_S <- if (grepl("^\\d+S", cigar)) as.integer(sub("^(\\d+)S.*", "\\1", cigar)) else 0L
    lead_H <- if (grepl("^\\d+H", cigar)) as.integer(sub("^(\\d+)H.*", "\\1", cigar)) else 0L
    lead   <- lead_S + lead_H
    
    # Trailing clips: soft (S) plus hard (H)
    trail_S <- if (grepl("\\d+S$", cigar)) as.integer(sub(".*?(\\d+)S$", "\\1", cigar)) else 0L
    trail_H <- if (grepl("\\d+H$", cigar)) as.integer(sub(".*?(\\d+)H$", "\\1", cigar)) else 0L
    trail   <- trail_S + trail_H
    
    c(lead, trail)
  }
  # *2) Helper: parse insertions and deletions from a CIGAR string
  .getIndels <- function(cig) {
    # Insertions in read‐space
    ins_ranges <- cigarRangesAlongQuerySpace(cig, ops="I")
    ins_pos    <- unlist(start(ins_ranges))
    ins_len    <- unlist(width(ins_ranges))
    
    # Deletions: pos in read‐space, len in reference‐space
    del_ranges_q <- cigarRangesAlongQuerySpace(cig, ops="D")
    del_ranges_r <- cigarRangesAlongReferenceSpace(cig, ops="D")
    del_pos      <- unlist(start(del_ranges_q))
    del_len      <- unlist(width(del_ranges_r))
    
    list(
      insertion = list(pos = ins_pos, len = ins_len),
      deletion  = list(pos = del_pos, len = del_len)
    )
  }
  
  # *3) Helper: find all cytosine positions in the read sequence
  .getCytPos <- function(read_seq, strand) {
    # If the read is reverse‐mapped, reverse‐complement it first
    seq_to_search <- if (strand == "-") reverseComplement(read_seq) else read_seq
    
    # Use matchPattern() to find every “C” and return their 1-based starts
    start(matchPattern("C", seq_to_search))
  }
  
  # Find the reference aligned insertion position using the cigar string 
  .getInsertionPos <- function(cigar){
    # Explode the CIGAR string into operations and lengths
    ops  <- explodeCigarOps(cigar)[[1]]
    lens <- explodeCigarOpLengths(cigar)[[1]]
    
    read_pos <- 1
    ref_pos  <- 1
    insertion_ref_sites <- c()
    
    for (i in seq_along(ops)) {
      op  <- ops[i]
      len <- lens[i]
      
      if (op == "M" || op == "=" || op == "X") {
        read_pos <- read_pos + len
        ref_pos  <- ref_pos + len
      } else if (op == "I") {
        insertion_ref_sites <- c(insertion_ref_sites, rep(ref_pos, len))
        read_pos <- read_pos + len
      } else if (op == "D" || op == "N") {
        ref_pos  <- ref_pos + len
      } else if (op == "S" || op == "H" || op == "P") {
        read_pos <- read_pos + len
      }
    }
    
    return(unique(insertion_ref_sites))
  }
  
  # *4) Main worker: decode MM/ML, apply soft‐clip and indel adjustments,
  #    and return a DataFrame of genomic positions + probabilities + modType
  
  .computeModifPos <- function(mm_tag, ml_tag, cigar, strand, start,
                               cyt_pos, indel_info, ref_seq, aligned_seq, 
                               context, synonymous) {
    
    # 0) bail out if no MM or ML tag
    if (is.na(mm_tag) || length(mm_tag)==0 ||
        is.null(ml_tag) || length(ml_tag)==0) {
      return(DataFrame(pos=integer(0),
                       prob=numeric(0),
                       modType=character(0)))
    }
    
    # 1) directly grab the ML vector of 0–255 probabilities
    ml_vals <- as.integer(ml_tag)
    ml_idx  <- 1L
    
    # 2) figure out read‐length and soft‐clipped bases
    read_len   <- cigarWidthAlongQuerySpace(cigar)
    read_ref_len <- cigarWidthAlongReferenceSpace(cigar) 
    sc         <- .getClipping(cigar)    # c(lead_clip, trail_clip)
    lead_clip  <- sc[1]
    trail_clip <- sc[2]
    
    # 3) pull out this read’s insertions/deletions
    ins         <- indel_info$insertion
    del         <- indel_info$deletion
    insert_pos  <- ins$pos;    insert_len  <- ins$len
    delete_pos  <- del$pos;    delete_len  <- del$len
    
    # filtering insertion position 
    ins_align_pos <- .getInsertionPos(cigar) - 1
    
    # 4) prepare outputs
    out_pos  <- integer()
    out_prob <- numeric()
    out_modified_base_code  <- character()
    out_mod <- integer()
    out_insert <- integer()
    out_delete <- integer()
    out_diff <- integer()
    out_nocall <-integer()
    
    
    # 5) split your “C+m?” / “C+h?” entries
    entries <- strsplit(as.character(mm_tag), ";", fixed=TRUE)[[1]]
    for (ent in entries) {
      if (ent == "") next  # skip empty
      
      parts    <- strsplit(ent, ",", fixed=TRUE)[[1]]
      code     <- parts[1]                     # e.g. "C+m?"
      offsets  <- as.integer(parts[-1]) + 1L   # 1‐based step counts
      c_counts <- cumsum(offsets)              # which nth‐C in the read
      
      # 6) grab the ML probabilities for *those* offsets
      probs    <- ml_vals[ ml_idx + seq_along(offsets) - 1L ] / 255
      ml_idx   <- ml_idx + length(offsets)
      
      # 7) throw away any c_counts beyond the # of Cs in the read
      valid <- c_counts <= length(cyt_pos)

      if (!any(valid)) next  
      
      # 8) apply your lead_clip window 
      hit_cyt_pos <- cyt_pos[c_counts[valid]]
      hit_probs <- probs[valid]
      hit_code <- rep(code, sum(valid))
      n_hits <- length(hit_cyt_pos)
      if (n_hits == 0) next
      
      # 9) for each one that survived, map read‐space→genome
      
      for (k in seq_along(hit_cyt_pos)) {
        # the read-space C position for this k
        pr <- hit_probs[k]
        rp <- hit_cyt_pos[k]
        cd <- hit_code[k]
        if (strand == "+") {
          # **collect any C that was part of an insertion**
          insert <- as.integer(any(rp >= insert_pos & rp < insert_pos + insert_len))
          # delete <- as.integer(any(rp == delete_pos))
          delete <- as.integer(any(rp %in% delete_pos))
          
          insert_off <- sum(insert_len[insert_pos < rp])
          delete_off <- sum(delete_len[delete_pos < rp])
          
          if (delete == 1){
            del_idx <- which(rp >= delete_pos & rp < delete_pos + 1)[1]
            del_len <- if (!is.na(del_idx)) delete_len[del_idx] else 0L
            gp <- start - lead_clip + (rp - insert_off + delete_off) - 1L + del_len
          } else {
            gp <- start - lead_clip + (rp - insert_off + delete_off) - 1L
          }
        } else {
          # reverse strand
          rp_rev <- read_len - rp + 1L
          
          insert <- as.integer(any(rp_rev >= insert_pos & rp_rev <= insert_pos + insert_len))
          # delete <- as.integer(any(rp_rev == delete_pos))
          delete <- as.integer(any(rp_rev %in% delete_pos))
          
          insert_off <- sum(insert_len[insert_pos < rp_rev])
          delete_off <- sum(delete_len[delete_pos < rp_rev])
          
          if (delete == 1) {
            del_idx <- which(rp_rev >= delete_pos & rp_rev < delete_pos + 1)[1]
            del_len <- if (!is.na(del_idx)) delete_len[del_idx] else 0L
            within_deletion <- any(rp_rev >= delete_pos & rp_rev < delete_pos + delete_len)
            adjust <- if (within_deletion) -1L else 0L
            gp <- start - lead_clip + (rp_rev - insert_off + delete_off) - 1L + del_len  + adjust
          } else {
            gp <- start - lead_clip + (rp_rev - insert_off + delete_off) - 1L
          }
        }
        if (gp >= start + read_ref_len | gp < start) next
        
        # Fetch ref/read base at this position
        ref_idx <- gp - start + 1
        if (ref_idx < 1 || ref_idx > length(ref_seq)) next
        ref_base <- ref_seq[ref_idx]
        aligned_base <- aligned_seq[ref_idx]
        
        # Diff flag: ref!=read (strict SNV filtering)
        if (strand == "+"){
          diff <- as.integer(ref_base != aligned_base | as.character(aligned_base) != "C")
        }else if(strand == "-"){
          diff <- as.integer(ref_base != aligned_base | as.character(aligned_base) != "G")
        }
        
        # filtering only if insertion is exactly *within* the C site (not adjacent)
        if (ref_idx %in% ins_align_pos & strand == "+") {
          insert = 1
        } else if ((ref_idx - 1) %in% ins_align_pos & strand == "-") {
          insert = 1
        }
        
        # Nocall flag: unmatched context (e.g. CG, CHG)
        nocall <- 0L
        if (context %in% "CG") {
          if (strand == "+") {
            ref_idx2 <- ref_idx + 1
            if (ref_idx2 >= 1 && ref_idx2 <= length(ref_seq)) {
              ref_next <- ref_seq[ref_idx2]
              aligned_next <- aligned_seq[ref_idx2]
              nocall <- as.integer(ref_next != aligned_next)
            }
          } else if (strand == "-") {
            ref_idx2 <- ref_idx - 1
            if (ref_idx2 >= 1 && ref_idx2 <= length(ref_seq)) {
              ref_prev <- ref_seq[ref_idx2]
              aligned_prev <- aligned_seq[ref_idx2]
              nocall <- as.integer(ref_prev != aligned_prev)
            }
          }
          
        } 
        if (context %in% "CHG") {
          if (strand == "+") {
            if (ref_idx + 2 <= length(ref_seq)) {
              ref_c  <- as.character(ref_seq[ref_idx])
              ref_h  <- as.character(ref_seq[ref_idx + 1])
              ref_g  <- as.character(ref_seq[ref_idx + 2])
              
              aln_c  <- as.character(aligned_seq[ref_idx])
              aln_h  <- as.character(aligned_seq[ref_idx + 1])
              aln_g  <- as.character(aligned_seq[ref_idx + 2])

              nocall <- 1L
              if (synonymous) {
                if (ref_c == "C" && aln_c == "C" &&
                    ref_g == "G" && aln_g == "G" &&
                    ref_h %in% c("A", "T", "C") && aln_h %in% c("A", "T", "C")) {
                  nocall <- 0L  # context preserved → allow
                }
              } else {
                # strict match
                if (ref_c == aln_c && ref_h == aln_h && ref_g == aln_g) {
                  nocall <- 0L
                }
              }
            }
          } else if (strand == "-") {
            if (ref_idx - 2 >= 1) {
              ref_g  <- as.character(ref_seq[ref_idx])
              ref_h  <- as.character(ref_seq[ref_idx - 1])
              ref_c  <- as.character(ref_seq[ref_idx - 2])
              
              aln_g  <- as.character(aligned_seq[ref_idx])
              aln_h  <- as.character(aligned_seq[ref_idx - 1])
              aln_c  <- as.character(aligned_seq[ref_idx - 2])
              
              nocall <- 1L
              if (synonymous) {
                if (ref_g == "G" && aln_g == "G" &&
                    ref_c == "C" && aln_c == "C" &&
                    ref_h %in% c("A", "T", "G") && aln_h %in% c("A", "T", "G")) {
                  nocall <- 0L
                }
              } else {
                if (ref_g == aln_g && ref_h == aln_h && ref_c == aln_c) {
                  nocall <- 0L
                }
              }
            }
          } 
        }
        if (context %in% "CHH"){
          if (strand == "+") {
            if (ref_idx + 2 <= length(ref_seq)) {
              ref_c  <- as.character(ref_seq[ref_idx])
              ref_h1  <- as.character(ref_seq[ref_idx + 1])
              ref_h2  <- as.character(ref_seq[ref_idx + 2])
              
              aln_c  <- as.character(aligned_seq[ref_idx])
              aln_h1  <- as.character(aligned_seq[ref_idx + 1])
              aln_h2  <- as.character(aligned_seq[ref_idx + 2])
              
              nocall <- 1L
              if (synonymous) {
                if (ref_c == "C" && aln_c == "C" &&
                    ref_h1 %in% c("A", "T", "C") && aln_h1 %in% c("A", "T", "C") &&
                    ref_h2 %in% c("A", "T", "C") && aln_h2 %in% c("A", "T", "C")) {
                  nocall <- 0L  # context preserved → allow
                }
              } else {
                # strict match
                if (ref_c == aln_c && ref_h1 == aln_h1 && ref_h2 == aln_h2) {
                  nocall <- 0L
                }
              }
            }
          } else if (strand == "-") {
            if (ref_idx - 2 >= 1) {
              ref_g  <- as.character(ref_seq[ref_idx])
              ref_h1  <- as.character(ref_seq[ref_idx - 1])
              ref_h2  <- as.character(ref_seq[ref_idx - 2])
              
              aln_g  <- as.character(aligned_seq[ref_idx])
              aln_h1  <- as.character(aligned_seq[ref_idx - 1])
              aln_h2  <- as.character(aligned_seq[ref_idx - 2])
              
              nocall <- 1L
              if (synonymous) {
                if (ref_g == "G" && aln_g == "G" &&
                    ref_h1 %in% c("A", "T", "G") && aln_h1 %in% c("A", "T", "G") &&
                    ref_h2 %in% c("A", "T", "G") && aln_h2 %in% c("A", "T", "G")) {
                  nocall <- 0L
                }
              } else {
                if (ref_g == aln_g && ref_h1 == aln_h1 && ref_h2 == aln_h2) {
                  nocall <- 0L
                }
              }
            }
          } 
        }
     
        # mod flag
        mod <- ifelse(insert == 0 & diff == 0 & nocall == 0, 1, 0)

        out_pos  <- c(out_pos,  gp)
        out_prob <- c(out_prob, pr)
        out_modified_base_code  <- c(out_modified_base_code, cd)
        out_mod <- c(out_mod, mod)
        out_insert <- c(out_insert, insert)
        out_delete <- c(out_delete, delete)
        out_diff <- c(out_diff, diff)
        out_nocall <- c(out_nocall, nocall)
        
        # # DEBUG
        # cat(sprintf("k=%d | start=%d | rp=%d | insert=%d | delete=%d |insert_off=%d | delete_off=%d | gp=%d | ref_base=%s | aln_base=%s\n  | del_len=%d",
        #             k, start, rp, insert, delete, insert_off, delete_off, gp, ref_base, aligned_base, del_len))
      }
      
    }
    DataFrame(pos=out_pos, prob=out_prob, modType=out_modified_base_code,
              Nmod=out_mod, Ninsert=out_insert, Ndelete=out_delete, Ndiff=out_diff, Nnocall=out_nocall)
  }
  
  # *5) Helper: Run .computeModifPos using the scanned bam list and return the modified
  # position and proportion values of the sequence as GRanges 
  .getCovSequence <- function(i) {
    if (i %in% progress_points) {
      pct <- round(i / N * 100)
      cat(sprintf("[readONTbam] Decoding MM/ML: %d%% (%d/%d reads)\n", pct, i, N))
    }
    # run .computeModifPos and generate the data_frame 
    df <- .computeModifPos(
      mm_tag      = sb$tag$MM[[i]],
      ml_tag      = sb$tag$ML[[i]],
      cigar       = sb$cigar[[i]],
      strand      = strands[[i]],
      start       = sb$pos[[i]],
      cyt_pos     = cyt_pos_list[[i]],
      indel_info  = indel_info_list[[i]],
      ref_seq     = ref_seq_list[[i]],
      aligned_seq = aligned_seq_list[[i]],
      context     = context,
      synonymous  = synonymous
    )
    
    # exclude the rows if they didn't hold the Nmod == 1
    df <- df[df$Nmod == 1, ]
    if (nrow(df) == 0) return(GRanges())
    
    # wrap up the results as the GRanges with using the index of sequence
    rr <- GRanges(
      seqnames = sb$rname[[i]],
      ranges   = IRanges(start = df$pos, width = 1),
      strand   = strands[[i]]
    )
    mcols(rr)$prob    <- df$prob
    mcols(rr)$modType <- df$modType
    mcols(rr)$read    <- read_order[[i]]
    
    # ## DEBUG
    # cat("current position i:", i, "\n")
    return(rr)
  }
  # -------------- 
  # *preparation of the decoding the modified information by
  # creating cytosine position and indel information list
  cat("[readONTbam] Generating cytosine position list and indel_info ...\n")
  cyt_pos_list <- mapply(.getCytPos,
                         read_seq = sb$seq,
                         strand   = strands,
                         SIMPLIFY = FALSE)
  
  indel_info_list   <- mapply(.getIndels,
                         cig      = sb$cigar,
                         SIMPLIFY = FALSE)
  
  # --------------    
  ## 2) decode every read
  cat("[readONTbam] Start Decoding MM/ML ...\n")
  N <- length(sb$qname)
  progress_points <- round(seq(0, N, length.out = 11)[-1])  # 10%, 20%, ..., 100%

  per_read <- bplapply(seq_len(N), .getCovSequence, BPPARAM = BPPARAM)
  
  valid_gr <- per_read[sapply(per_read, function(x) is(x, "GRanges") && length(x) > 0)]
  
  if (length(valid_gr) == 0) {
    all_calls <- GRanges()
  } else {
    all_calls <- suppressWarnings(do.call(c, valid_gr))
  }
  rm(per_read)
  
  
  # @@@@ time t4
  t4 <- proc.time()
  cat(sprintf("Step 3 (decode Bam file): %.2f sec\n", (t4-t3)[["elapsed"]]))
  
  ## 3) initialize ONT_Cm and ONT_C as simple character columns
  mcols(ref_gr)$ONT_Cm <- CharacterList(vector("list", length(ref_gr)))
  mcols(ref_gr)$ONT_C  <- CharacterList(vector("list", length(ref_gr)))
  
  ref_gr$readsM <- integer(length(ref_gr))
  ref_gr$readsN <- integer(length(ref_gr))
  
  
  cat("[readONTbam] Calculating overlaps with reference cytosines ...\n")
  ov_cov <- findOverlaps(ref_gr, all_calls, ignore.strand = FALSE)
  df_cov <- DataFrame(
    qhit   = queryHits(ov_cov),
    read   = all_calls$read[subjectHits(ov_cov)],
    rstrand= as.character(strand(all_calls))[subjectHits(ov_cov)],
    modType = mcols(all_calls)$modType[subjectHits(ov_cov)],
    prob    = mcols(all_calls)$prob[subjectHits(ov_cov)]
  )
  
  # subsetting the df_cov as the modif preset
  df_cov   <- df_cov[df_cov$modType %in% modif, ]
  
  # @@@@ time t5
  t5 <- proc.time()
  cat(sprintf("Step 4 (coverage): %.2f sec\n", (t5-t4)[["elapsed"]]))
  
  ## 4)  save the modified and read‐indices of modified (readsN, ONT_C) in ref_gr
  cat("[readONTbam] Annotating high-confidence modified reads ...\n")

  df_mod   <- df_cov[df_cov$modType %in% modif & df_cov$prob >= prob_thresh, ]
  mod_list <- lapply(split(df_mod$read, df_mod$qhit), unique)
  mod_char_list <- lapply(mod_list, function(x) paste(x, collapse = ","))
  # collapse into semicolon-separated strings
  idx_mod <- as.integer(names(mod_list))
  
  ref_gr$ONT_Cm[idx_mod] <- CharacterList(mod_char_list)
  ref_gr$readsM[idx_mod] <- lengths(mod_list)
  
  # @@@@ time t6
  t6 <- proc.time()
  cat(sprintf("Step 5 (ONT_Cm & readsM meta-column): %.2f sec\n", (t6-t5)[["elapsed"]]))
  
  # 5) save the coverage and read‐indices of unmodified (readsN, ONT_C) in ref_gr
  cat("[readONTbam] Annotating unmodified (covered) reads & Calculating coverage per position...\n")  
  
  df_unmod   <- df_cov[df_cov$modType %in% modif & df_cov$prob < prob_thresh, ]
  
  unmod_list <- lapply(split(df_unmod$read, df_unmod$qhit), unique)
  cov_list   <- lapply(split(df_cov$read, df_cov$qhit), unique)
  unmod_char_list <- lapply(unmod_list, function(x) paste(x, collapse = ","))
  # collapse into semicolon-separated strings
  idx_unmod <- as.integer(names(unmod_list))
  idx_cov   <- as.integer(names(cov_list))
  
  ref_gr$ONT_C[idx_unmod] <- CharacterList(unmod_char_list)
  ref_gr$readsN[idx_cov] <- lengths(cov_list)
  
  # @@@@ time t6
  t7 <- proc.time()
  cat(sprintf("Step 6 (ONT_C & readsN meta-column): %.2f sec\n", (t7-t6)[["elapsed"]]))
  cat("[readONTbam] Done. Returning annotated GRanges.\n")
  cat(sprintf("Total: %.2f sec\n", (t7-t1)[["elapsed"]]))
  
  return(ref_gr)
  
}

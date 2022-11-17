#'
#' Sorted Annotations
#'
#' @description Sorts annotations by position and chromosome.
#'
#' @param granges An object of GRanges or data.frame type. For more details
#' please check ?GRanges.sortAnnotations and ?data.table.sortAnnotations.
#' @param chromosomes A list of chromosomes which should be returned from the
#' full annotations. Returns all chromosomes if NULL.
#' @param verbose Should workflow messages be printed? (default: FALSE)
#' @param logFunction Function to use for printing messages. Useful if you
#' wish to setup custom logging (default: message).
#'
#' @importFrom data.table as.data.table
#'
sortAnnotations <- function(granges, chromosomes = NULL, verbose = FALSE, logFunction = message) {
  if ('GRanges' %in% class(granges)) {
    return(GRanges.sortAnnotations(granges, chromosomes = chromosomes, verbose = verbose, logFunction = logFunction))
  }

  if (any(c('data.frame', 'data.table', 'DFrame') %in% class(granges))) {
    return(data.table.sortAnnotations(data.table::as.data.table(granges), chromosomes = chromosomes, verbose = verbose, logFunction = logFunction))
  }

  stop('Annotation sorting not implemented for objects of class "', paste(class(granges), collapse = '", "'), , '". Please convert it to GRanges, data.frame or data.table object.')
}

#'
#' Sorted Annotations
#'
#' @description Sorts annotations by position and chromosome.
#'
#' @param granges A GRanges object. Should contain CpG IDs as names or have a
#' "Name" metadata column.
#' @param chromosomes A list of chromosomes which should be returned from the
#' full annotations. Returns all chromosomes if NULL.
#' @param verbose Should workflow messages be printed? (default: FALSE)
#' @param logFunction Function to use for printing messages. Useful if you
#' wish to setup custom logging (default: message).
#'
#' @importFrom GenomicRanges seqnames sort mcols
#'
GRanges.sortAnnotations <- function(granges, chromosomes = NULL, verbose = FALSE, logFunction = message) {
  if (is.null(chromosomes)) {
    if (verbose) logFunction('[ Selecting all available chromosomes from the supplied annotation ]')
    chromosomes <- GenomicRanges::seqnames(granges)
  } else {
    if (verbose) logFunction(paste0('[ Selecting ', paste(chromosomes, collapse = ', '), ' chromosomes from the supplied annotation ]'))
  }

  if (is.null(names(granges))) {
    if (!('Name' %in% colnames(GenomicRanges::mcols(granges)))) {
      stop('Unable to determine probe positions from current annotation. Please set GRanges names to probe names or add a "Name" metadata column to your GRanges object')
    }

    if (verbose) logFunction('[ Setting GRanges names to probe names ]')
    names(granges) <- GenomicRanges::mcols(granges)[ , 'Name' ]
  }

  granges <- granges[ GenomicRanges::seqnames(granges) %in% chromosomes ]
  granges <- GenomicRanges::sort(granges)

  return(granges)
}

#'
#' Sorted Annotations
#'
#' @description Sorts annotations by position and chromosome.
#'
#' @param granges A data.table object containing CpG annotations. Should contain
#' 'Name', 'chr' and 'pos' columns.
#' @param chromosomes A list of chromosomes which should be returned from the
#' full annotations. Returns all chromosomes if NULL.
#' @param verbose Should workflow messages be printed? (default: FALSE)
#' @param logFunction Function to use for printing messages. Useful if you
#' wish to setup custom logging (default: message).
#'
#' @importFrom data.table .I
#' @importFrom dplyr %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols sort
#'
data.table.sortAnnotations <- function(granges, chromosomes = NULL, verbose = FALSE, logFunction = message) {
  if (!all(c('Name', 'chr', 'pos') %in% colnames(granges))) {
    stop('Missing columns in annotation. These columns should be present: "Name", "chr" and "pos"')
  }

  if (is.null(chromosomes)) {
    if (verbose) logFunction('[ Selecting all available chromosomes from the supplied annotation ]')
    chromosomes <- unique(granges[ , chr ])
  } else {
    if (verbose) logFunction(paste0('[ Selecting ', paste(chromosomes, collapse = ', '), ' chromosomes from the supplied annotation ]'))
  }

  map <- granges %>%
    .[ , list(chr, pos, Name) ] %>%
    .[ chr %in% chromosomes ] %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, start.field = 'pos', end.field = 'pos')

  # Update rownames of the granges object
  names(map) <- GenomicRanges::mcols(map) %>% .[ , 1 ]
  GenomicRanges::mcols(map) <- NULL
  map <- GenomicRanges::sort(map)

  return(map)
}

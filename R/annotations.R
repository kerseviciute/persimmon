#'
#' Sorted annotations
#'
#' @description Get sorted annotations from _RGChannelSet_, _RGChannelSetExtended_ or _MethylSet_
#' object from specified chromosomes.
#'
#' @param set a _RGChannelSet_, _RGChannelSetExtended_ or _MethylSet_ object
#' @param chromosomes a list of chromosomes for which to extract annotations. Will
#' return all chromosomes if not specified.
#'
#' @import minfi
#' @import data.table
#' @import GenomicRanges
#'
#' @export
#'
extractSortedAnnotations <- function(set, chromosomes = NULL) {
  classes <- c('RGChannelSet', 'RGChannelSetExtended', 'MethylSet')
  if (!any(classes %in% class(set))) {
    stop('Object should be of type ', paste(classes, collapse = ', '))
  }

  annotations <- minfi::getAnnotation(set) %>%
    as.data.table

  if (is.null(chromosomes)) {
    chromosomes <- unique(annotations[ , chr ])
  }

  map <- annotations %>%
    .[ , list(chr, pos, Name) ] %>%
    .[ chr %in% chromosomes ] %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, start.field = 'pos', end.field = 'pos')

  # Update rownames of the granges object
  names(map) <- GenomicRanges::mcols(map) %>% .[ , 1 ]
  GenomicRanges::mcols(map) <- NULL
  map <- GenomicRanges::sort(map)

  return(map)
}

#'
#' @export
#'
sortAnnotations <- function(granges, ...) {
  if ('GRanges' %in% class(granges)) {
    return(sortAnnotationsGRanges(granges, ...))
  }

  if (any(c('data.frame', 'data.table', 'DFrame') %in% class(granges))) {
    return(sortAnnotationsDataTable(as.data.table(granges), ...))
  }

  stop('Annotation sorting not implemented for objects of class "', paste(class(granges), collapse = '", "'), , '". Please convert it to GRanges, data.frame or data.table object.')
}

#'
#' @export
#'
sortAnnotationsGRanges <- function(granges, chromosomes = NULL, verbose = FALSE) {
  if (is.null(chromosomes)) {
    if (verbose) message('[ Selecting all available chromosomes from the supplied annotation ]')
    chromosomes <- seqnames(granges)
  } else {
    if (verbose) message('[ Selecting ', paste(chromosomes, collapse = ', '), ' chromosomes from the supplied annotation ]')
  }

  if (is.null(names(granges))) {
    if (!('Name' %in% colnames(mcols(granges)))) {
      stop('Unable to determine probe positions from current annotation. Please set GRanges names to probe names or add a "Name" metadata column to your GRanges object')
    }

    if (verbose) message('[ Setting GRanges names to probe names ]')
    names(granges) <- mcols(granges)[ , 'Name' ]
  }

  granges <- granges[ seqnames(granges) %in% chromosomes ]
  granges <- GenomicRanges::sort(granges)

  return(granges)
}

#'
#' @export
#'
sortAnnotationsDataTable <- function(granges, chromosomes = NULL, verbose = FALSE) {
  if (!all(c('Name', 'chr', 'pos') %in% colnames(granges))) {
    # TODO: proper error message
    stop('Missing columns')
  }

  if (is.null(chromosomes)) {
    if (verbose) message('[ Selecting all available chromosomes from the supplied annotation ]')
    chromosomes <- unique(granges[ , chr ])
  } else {
    if (verbose) message('[ Selecting ', paste(chromosomes, collapse = ', '), ' chromosomes from the supplied annotation ]')
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

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



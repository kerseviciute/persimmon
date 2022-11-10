#'
#' Sorted annotations
#'
#' @description Get sorted annotations from _RGChannelSet_ or _RGChannelSetExtended_
#' object from specified chromosomes.
#'
#' @param rgset a _RGChannelSet_ or _RGChannelSetExtended_ object
#' @param chromosomes a list of chromosomes for which to extract annotations. Will
#' return all chromosomes if not specified.
#'
#' @import minfi
#' @import data.table
#' @import GenomicRanges
#'
#' @export
#'
extractSortedAnnotations <- function(rgset, chromosomes = NULL) {
  if (!any(class(rgset) %in% c('RGChannelSet', 'RGChannelSetExtended'))) {
    stop('Expecting an object of type \'RGChannelSet\' or \'RGChannelSetExtended\'')
  }

  annotations <- minfi::getAnnotation(rgset) %>%
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



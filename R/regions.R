#'
#' GRanges Splitting
#'
#' @description Split GRanges of a single chromosome into splits smaller than maxSplit size.
#' Takes into account the distance between individual CpGs and performs splitting at
#' locations which are the furthest from each other.
#'
#' @param granges GRanges object
#' @param maxSplit maximum split size, integer
#'
#' @import GenomicRanges
#'
split <- function(granges, maxSplit = 50000) {
  if (length(granges) < maxSplit) {
    return(granges)
  }

  delta <- diff(start(granges))
  splitAt <- which.max(delta)

  grangesA <- granges[ 1:splitAt ]
  grangesB <- granges[ (splitAt + 1):length(granges) ]

  if (length(grangesA) > maxSplit) {
    grangesA <- split(grangesA)
  }

  if (length(grangesB) > maxSplit) {
    grangesB <- split(grangesB)
  }

  list(grangesA, grangesB)
}

#'
#' GRanges Splitting
#'
#' @description Split GRanges object into splits smaller than maxSplit size.
#' Each chromosome is split separately, at locations which are furthest apart.
#'
#' @param granges GRanges object
#' @param maxSplit maximum split size, integer
#'
#' @import GenomicRanges
#'
splitChromosomes <- function(granges, maxSplit = 50000) {
  byChromosome <- GenomicRanges::split(granges, seqnames(granges))

  foreach(chr = byChromosome) %do% {
    split(chr, maxSplit = maxSplit)
  } %>%
    unlist %>%
    GenomicRanges::GRangesList()
}

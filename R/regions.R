#'
#' GRanges Splitting
#'
#' @description Splits GRanges of a single chromosome into splits smaller than maxSplit size.
#' Takes into account the distance between individual CpGs and splits the dataset
#' between CpGs that are furthest apart.
#'
#' @param granges GRanges object.
#' @param maxSplit Maximum split size, integer (default: 50000).
#'
#' @importFrom GenomicRanges start
#'
split <- function(granges, maxSplit = 50000) {
  if (length(granges) < maxSplit) {
    return(granges)
  }

  delta <- diff(GenomicRanges::start(granges))
  splitAt <- which.max(delta)

  grangesA <- granges[ 1:splitAt ]
  grangesB <- granges[ (splitAt + 1):length(granges) ]

  if (length(grangesA) > maxSplit) {
    grangesA <- split(grangesA, maxSplit)
  }

  if (length(grangesB) > maxSplit) {
    grangesB <- split(grangesB, maxSplit)
  }

  list(grangesA, grangesB)
}

#'
#' GRanges Splitting
#'
#' @description Splits GRanges object into splits smaller than maxSplit size.
#' Each chromosome is split separately, between CpGs that are furthest apart.
#'
#' @param granges GRanges object.
#' @param maxSplit Maximum split size, integer (default: 50000).
#'
#' @importFrom GenomicRanges split seqnames GRangesList
#' @importFrom foreach foreach %do%
#' @importFrom dplyr %>%
#'
splitChromosomes <- function(granges, maxSplit = 50000) {
  byChromosome <- GenomicRanges::split(granges, GenomicRanges::seqnames(granges))

  splits <- foreach::foreach(chr = byChromosome) %do% {
    split(chr, maxSplit = maxSplit)
  } %>%
    unlist %>%
    GenomicRanges::GRangesList()

  names(splits) <- make.names(seq_along(splits))

  splits
}

#'
#' Methylation Regions
#'
#' @description Calculates methylation regions for a subset of the dataset.
#' When using this method, a single chromosome or its subset should be passed.
#'
#' @return Returns a named list. The names correspond to CpG names and the
#' values are for that CpG assigned cluster. Note that not all CpGs that
#' were in the original dataset will be returned. Some may not cluster at
#' all, and thus will not be returned.
#'
#' @param beta Beta value (or m-value) matrix with probes in rows and
#' samples in columns. Row names should correspond to names of the annotation
#' (see \code{granges} argument).
#' @param granges CpG annotations. It is a GRanges object containing sorted
#' coordinates of all CpGs of interest. The names of this object should
#' correspond to CpG names and match row names of the beta value matrix.
#' @param thresholdDistance Threshold distance
#' @param rateFactor Rate factor
#' @param scale Should beta value matrix be scaled? (default: TRUE)
#' @param minClusterSize Minimum number of CpGs per cluster, integer
#' (default: 5).
#' @param compressionRatio Compression ratio
#' @param seed A seed for random number generator.
#' @param verbose Should workflow messages be printed? (default: FALSE)
#'
#' @importFrom Rfast Dist
#' @importFrom stats ecdf
#' @importFrom stats hclust cutree
#'
methRegions <- function(
  beta,
  granges,
  thresholdDistance = 5e+07,
  rateFactor = 10,
  scale = TRUE,
  minClusterSize = 5,
  compressionRatio = 10,
  seed = 23985034,
  verbose = FALSE) {
  set.seed(seed)

  if (!all(names(granges) %in% rownames(beta))) {
    warning('Some probes are missing in the beta value beta matrix. Will use only positions that are present in both annotation and beta value matrix.')

    beta <- beta[ na.omit(match(names(granges), rownames(beta))), ]
    granges <- granges[ names(granges) %in% rownames(beta) ]
  }

  beta <- beta[ names(granges), ]
  if (anyNA(beta)) {
    stop('Methylation values should not contain missing values')
  }

  # TODO: discuss when this step should be performed
  if (verbose) message('[ Preparing methylation matrix ]')
  beta <- prepareBeta(beta, scale = scale, verbose = verbose)

  if (verbose) message('[ Calculating euclidean distance between probe methylation levels ]')
  distance <- Rfast::Dist(beta)
  distance <- distance / max(distance, na.rm = TRUE)
  colnames(distance) <- rownames(distance) <- names(granges)

  # TODO: this needs to be checked
  w <- base::mapply(function(x) { x - start(granges) }, start(granges))
  w <- abs(w)
  colnames(w) <- rownames(w) <- names(granges)

  if (verbose) message('[ Estimating sigmoid function center ]')
  n <- min(1e7, length(w))
  percentile <- stats::ecdf(sample(w, n))(thresholdDistance)
  center <- percentile
  rate <- 1 / percentile * rateFactor

  if (verbose) message('[ Computing scaled distance matrix ]')
  ww <- w / max(w)
  ww <- sigmoid(ww - center, k = rate)
  X <- (1 - ww) * distance + (ww) * 1
  diag(X) <- 0

  if (verbose) message('[ Clustering scaled distance matrix ]')
  hclustering <- stats::hclust(as.dist(X), method = 'complete')

  if (verbose) message('[ Cutting tree and calculating clusters ]')
  # TODO: test dynamic tree cut
  k <- max(floor(ncol(X) / compressionRatio), 1)
  clusters <- stats::cutree(hclustering, k = k)
  bad <- which(table(clusters) < minClusterSize)
  clusters[ clusters %in% bad ] <- 0
  names(clusters) <- names(granges)
  clusters <- clusters[ clusters != 0 ]

  return(clusters)
}

#'
#' Sigmoid function
#'
#' @param x Matrix
#' @param k k
#'
sigmoid <- function(x, k = 1) {
  1 / (1 + exp(-k * x))
}

#'
#' Beta Values
#'
#' @description Beta value preprocessing.
#'
#' @param beta Beta value matrix.
#' @param scale Should beta value matrix be scaled? (default: TRUE)
#' @param verbose Should workflow messages be printed? (default: FALSE)
#'
prepareBeta <- function(beta, scale = TRUE, verbose = TRUE) {
  if (scale == TRUE) {
    if (verbose) message('[ Scaling methylation matrix ]')
    beta <- beta %>% base::t() %>% base::scale() %>% base::t()
  } else {
    if (verbose) message('[ Using non-scaled methylation matrix ]')
  }

  return(beta)
}

#'
#' Methylation Regions
#'
#' @inheritParams methRegions
#' @inheritParams sortAnnotations
#' @inheritParams splitChromosomes
#' @param allowParallel Should parallel processing be used. If you wish
#' to use multiple cores, register cores using \code{doParallel} package.
#'
#' @importFrom tibble enframe
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom dplyr %>%
#' @importFrom data.table .I as.data.table ':=' .N data.table setnames
#'
#' @export
#'
findMethRegions <- function(
  beta,
  granges,
  chromosomes = NULL,
  maxSplit = 50000,
  thresholdDistance = 5e+07,
  rateFactor = 10,
  scale = TRUE,
  minClusterSize = 5,
  compressionRatio = 10,
  seed = 23985034,
  allowParallel = FALSE,
  verbose = FALSE
) {
  granges <- sortAnnotations(granges, verbose = verbose, chromosomes = chromosomes)
  granges <- splitChromosomes(granges, maxSplit = maxSplit)
  if (verbose) message('[ Dataset was divided into ', length(granges), ' splits ]')

  if (allowParallel) {
    if (verbose) message('[ Will try to parallelize using doParallel ]')

    clusters <- foreach::foreach(i = names(granges), .combine = rbind) %dopar% {
      methRegions(beta, granges[[ i ]], verbose = verbose,
                  thresholdDistance = thresholdDistance,
                  rateFactor = rateFactor, scale = scale,
                  minClusterSize = minClusterSize,
                  compressionRatio = compressionRatio, seed = seed) %>%
        tibble::enframe(name = 'Name', 'Cluster') %>%
        data.table::as.data.table %>%
        .[ , Cluster := paste(i, Cluster, sep = '_') ]
    }
  } else {
    if (verbose) message('[ Parallel processing is disabled ]')

    clusters <- foreach::foreach(i = names(granges), .combine = rbind) %do% {
      methRegions(beta, granges[[ i ]], verbose = verbose,
                  thresholdDistance = thresholdDistance,
                  rateFactor = rateFactor, scale = scale,
                  minClusterSize = minClusterSize,
                  compressionRatio = compressionRatio, seed = seed) %>%
        tibble::enframe(name = 'Name', 'Cluster') %>%
        data.table::as.data.table %>%
        .[ , Cluster := paste(i, Cluster, sep = '_') ]
    }
  }

  # Correct cluster names
  clusters <- data.table::data.table(Old = clusters[ , unique(Cluster) ]) %>%
    .[ , New := make.names(1:.N) ] %>%
    merge(clusters, by.x = 'Old', by.y = 'Cluster') %>%
    data.table::setnames('New', 'ClusterID') %>%
    .[ , Old := NULL ] %>%
    .[ , list(Name, ClusterID) ]

  clusters
}

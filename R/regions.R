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
#' TODO: explain all variables
#' TODO: create an explanatory page with examples how to run the package
#' TODO: check all calculations
#' TODO: is there a way to optimize?
#' TODO: what could go wrong? at what scenarios will the sigmoid fail?
#' TODO: should this method be exported?
#' TODO: can it be used with m-values, not beta values? why? explanations!
#'
#' @importFrom Rfast Dist
#' @importFrom stats ecdf
#' @importFrom stats hclust
#'
#'
#' @export
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
    notPresent <- names(granges)[ !(names(granges) %in% rownames(beta)) ]
    stop('Could not find probes ', paste(notPresent, collapse = ', '), ' in the beta value beta matrix')
  }

  beta <- beta[ names(granges), ]
  if (anyNA(beta)) {
    stop('Methylation values should not contain missing values')
  }

  # TODO: move out to wrapper
  if (verbose) message('[ Preparing methylation matrix ]')
  beta <- prepareBeta(beta, scale = scale, verbose = verbose)

  if (verbose) message('[ Calculating euclidean distance between probe methylation levels ]')
  distance <- Rfast::Dist(beta)
  distance <- distance / max(distance, na.rm = TRUE)
  colnames(distance) <- rownames(distance) <- names(granges)

  # TODO: this needs to be checked
  w <- mapply(function(x) { x - start(granges) }, start(granges))
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
  clusters <- cutree(hclustering, k = k)
  bad <- which(table(clusters) < minClusterSize)
  clusters[ clusters %in% bad ] <- 0
  names(clusters) <- names(granges)
  clusters <- clusters[ clusters != 0 ]

  return(clusters)
}

#'
#' @export
#'
sigmoid <- function(x, k = 1) {
  1 / (1 + exp(-k * x))
}

#'
#' @export
#'
prepareBeta <- function(beta, scale = TRUE, verbose = TRUE) {
  if (scale == TRUE) {
    if (verbose) message('[ Scaling methylation matrix ]')
    beta <- beta %>% t %>% scale %>% t
  } else {
    if (verbose) message('[ Using non-scaled methylation matrix ]')
  }

  return(beta)
}

#'
#' @importFrom tibble enframe
#' @import foreach
#' @import dplyr
#' @import data.table
#'
#'
#' @export
#'
findMethRegions <- function(granges, beta, allowParallel = FALSE, verbose = FALSE, chromosomes = NULL, maxSplit = 50000, ...) {
  granges <- sortAnnotations(granges, verbose = verbose, chromosomes = chromosomes)
  granges <- splitChromosomes(granges, maxSplit = maxSplit)
  if (verbose) message('[ Dataset was divided into ', length(granges), ' splits ]')

  if (allowParallel) {
    if (verbose) message('[ Will try to parallelize using doParallel ]')

    clusters <- foreach::foreach(i = names(granges), .combine = rbind) %dopar% {
      methRegions(beta, granges[[ i ]], verbose = verbose, ...) %>%
        tibble::enframe(name = 'Name', 'Cluster') %>%
        as.data.table %>%
        .[ , Cluster := paste(i, Cluster, sep = '_') ]
    }
  } else {
    if (verbose) message('[ Parallel processing is disabled ]')

    clusters <- foreach::foreach(i = names(granges), .combine = rbind) %do% {
      methRegions(beta, granges[[ i ]], verbose = verbose, ...) %>%
        tibble::enframe(name = 'Name', 'Cluster') %>%
        as.data.table %>%
        .[ , Cluster := paste(i, Cluster, sep = '_') ]
    }
  }

  # Correct cluster names
  clusters <- data.table(Old = clusters[ , unique(Cluster) ]) %>%
    .[ , New := make.names(1:.N) ] %>%
    merge(clusters, by.x = 'Old', by.y = 'Cluster') %>%
    setnames('New', 'ClusterID') %>%
    .[ , Old := NULL ] %>%
    .[ , list(Name, ClusterID) ]

  clusters
}

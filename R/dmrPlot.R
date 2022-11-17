#'
#' DMR plot
#'
#' @description Plots methylation region. Displays beta values for each same and
#' draws a curve across the mean group methylation value of each position.
#'
#' @return A \code{ggplot2} object.
#'
#' @param beta A beta value matrix.
#' @param cpgs A list of CpGs to visualize.
#' @param annotation CpG annotations containing position coordinates and CpG names.
#' @param conditions A named list of conditions. The values should contain sample groups
#' while the names should be sample IDs (same as colnames of the beta value matrix).
#'
#' @importFrom ggalt geom_xspline
#' @importFrom ggplot2 ggplot geom_point aes theme element_text
#' @importFrom data.table .I .BY ':=' setDT as.data.table
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @importFrom tibble enframe
#'
#' @export
#'
plotDMR <- function(beta, cpgs, annotation, conditions) {
  stopifnot(all(c('Name', 'pos') %in% colnames(annotation)))
  stopifnot(all(names(conditions) %in% colnames(beta)))
  stopifnot(all(cpgs %in% rownames(beta)))
  stopifnot(all(cpgs %in% annotation[ , Name ]))

  key <- conditions %>%
    tibble::enframe(name = 'SID', value = 'Condition')

  beta <- beta[ cpgs, ] %>%
    t %>%
    as.data.table(keep.rownames = 'SID') %>%
    merge(key, by = 'SID') %>%
    reshape2::melt(id.vars = c('SID', 'Condition'), variable.name = 'Name', value.name = 'Beta') %>%
    setDT %>%
    na.omit %>%
    merge(annotation[ , list(Name, Position = pos) ], by = 'Name') %>%
    .[ , Position := factor(Position, levels = sort(unique(Position), decreasing = FALSE)) ] %>%
    .[ order(Position) ]

  meanBeta <- beta %>%
    .[ , list(Beta = mean(Beta)), by = list(Condition, Position, Name) ] %>%
    .[ order(Position, decreasing = FALSE) ]

  ggplot2::ggplot() +
    ggplot2::geom_point(data = beta, ggplot2::aes(x = Position, y = Beta, color = Condition), alpha = 0.5) +
    ggalt::geom_xspline(data = meanBeta, ggplot2::aes(x = Position, y = Beta, group = Condition, color = Condition), spline_shape = 0.3, open = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

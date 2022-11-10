#'
#' Extract Beta or M-values
#'
#' @description Get beta or mvalues from _RGChannelSet_, _RGChannelSetExtended_ or _MethylSet_ object.
#'
#' @param rgset a _RGChannelSet_, _RGChannelSetExtended_ or _MethylSet_ object.
#' @param type what type of values should be used (_beta_ or _mvalues_)
#'
#' @import minfi
#'
getValues <- function(set, type = c('beta', 'mvalues')) {
  classes <- c('RGChannelSet', 'RGChannelSetExtended', 'MethylSet')
  if (!any(classes %in% class(set))) {
    stop('Object should be of type ', paste(classes, collapse = ', '))
  }

  type <- match.arg(type)

  if (class(set) %in% c('RGChannelSet', 'RGChannelSetExtended') & type == 'mvalues') {
    set <- minfi::preprocessRaw(set)
  }

  values <- switch(type,
                   'beta' = minfi::getBeta(set),
                   'mvalues' = minfi::getM(set))

  return(values)
}

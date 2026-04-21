
#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is to generate default knots
#' This code belongs to AdaptFitOS
#' @keywords internal
#'
#'

arg.searchOS <- function (string, arg.name)
{
  out <- SemiPar::break.string(string, arg.name)
  if (length(out) == 1)
    present <- FALSE
  if (length(out) > 1)
    present <- TRUE
  if (present) {
    right.string <- out[2]
    type.arg <- "ordinary"
    comma.found <- FALSE
    for (i in 1:nchar(right.string)) {
      if (substring(right.string, i, i) == ",")
        comma.found <- TRUE
      if ((substring(right.string, i, i) == "(") & comma.found ==
          FALSE)
        type.arg <- "array"
    }
    if (type.arg == "ordinary") {
      out.comma <- SemiPar::break.string(right.string, ",")
      arg.assign <- paste(arg.name, out.comma[1], sep = "")
    }
    if (type.arg == "array") {
      out.left <- SemiPar::break.string(right.string, "(")
      out <- SemiPar::break.string(right.string, ")")
      arg.assign <- paste(arg.name, out[1], ")", sep = "")
    }
  }
  if (!present)
    arg.assign <- NULL
  return(list(arg = arg.assign, present = present))
}

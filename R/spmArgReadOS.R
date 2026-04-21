#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is to generate default knots
#' This code belongs to AdaptFitOS
#' @keywords internal
#'
#'

spmArgReadOS <- function (arg.assignment)
{
    out <- SemiPar::break.string(arg.assignment, "=")
    arg.name <- out[1]
    arg.val <- eval(parse(text = out[2]), envir = sys.frame(-4))
    return(list(name = arg.name, val = arg.val))
}

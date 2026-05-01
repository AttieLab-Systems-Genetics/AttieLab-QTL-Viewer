#' Null-coalescing operator
#'
#' Returns left-hand side if not NULL, otherwise right-hand side.
#' @param a any
#' @param b any
#' @return a if not NULL else b
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

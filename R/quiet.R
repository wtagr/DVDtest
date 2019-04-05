#' Quite functions
#' 
#' @param x the stuff need to be quiet, messages of which are removed
#' @return quite result
#' @author Meng Xu
#' @keywords internal
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
### A collection of functions I find useful in my work as a genomic scientist.


#' Counts per million normalization
#'
#' This function converts an unnormalized countmatrix of genes x samples into counts-per-million normalization
#' @param data The count matrix to be normalized
#' @keywords normalization
#' @export
#' @examples
#' cpm()

cpm = function(data)
{
  data = apply(data, 2, function(x) (x/sum(x))*10^6)
}

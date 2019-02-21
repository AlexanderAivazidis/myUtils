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
  return(data)
}

#' Counts per million normalization
#'
#' This function converts an unnormalized countmatrix of genes x samples into read-per-kilobase normalization
#' @param data The count matrix to be normalized
#' @param geneLengths The gene lengths in the samer order as the rows of the count matrix
#' @keywords normalization
#' @export
#' @examples
#' rpkm()

rpkm = function(data,geneLengths){
  totalNumReads = colSums(data)
  test = data / sapply(1:length(totalNumReads), function(x) geneLengths/10^3 * totalNumReads[x]/10^6)
  return(data)
}





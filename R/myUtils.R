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
  data = data / sapply(1:length(totalNumReads), function(x) geneLengths/10^3 * totalNumReads[x]/10^6)
  return(data)
}

#' Map gene ids between one another for mouse
#' 
#' This function maps ids such as symbols , entrez ids, ensembls ids between each other
#'
#' @param IDs The IDs to be mapped
#' @param IDFrom The current ID format, such as SYMBOL, ENSEMBL
#' @param IDTo The ID type to map to, such as SYMBOL, ENSEMBL
#' @export

mapIdsMouse<-function(IDs,IDFrom,IDTo){
  require(org.Mm.eg.db)
  idmap=mapIds(x = org.Mm.eg.db,keys = IDs,column = IDTo, keytype = IDFrom,multiVals = "first")
  na_vec=names(idmap[is.na(idmap)==T])
  idmap=idmap[is.na(idmap)==F]
  idmap_df=data.frame("From"=names(idmap),"To"=unlist(unname(idmap)),stringsAsFactors = F)
  idmap_df = idmap_df[match(idmap_df[,2], idmap_df[,2]),]
  idmap_df = idmap_df[match(idmap_df[,1], idmap_df[,1]),]
  IDs[match(idmap_df[,1], IDs)] = idmap_df[,2]
  return(IDs)
}




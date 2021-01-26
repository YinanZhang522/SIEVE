#' Find the high frequency genes of gene sets defined by different method.
#' @name SIEVE
#' @usage SIEVE(genelist,n=50)
#' @param genelist The gene sets matrix. Rows should be genes and columns should be different batches.
#' @param n The minimum number of repetitions of a gene.


#' @return A list of high frequency genes.
#' @export

#' @examples geneset<-SIEVE(genelist,n=50)
SIEVE<-function(genelist,n=50){
  a_res<-table(as.vector(genelist))
  geneset<-names(a_res[a_res==n])
  return(geneset)
}

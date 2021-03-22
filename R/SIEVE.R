#' Find the high frequency genes of gene sets defined by different method.
#' @name SIEVE
#' @usage SIEVE(genelist,n=50)
#' @param genelist The gene sets matrix. Rows should be genes and columns should be different batches.
#' @param n The minimum number of repetitions of a gene.


#' @return A list of high frequency genes and a schematic diagram of the gene frequency.
#' @export

#' @examples geneset<-SIEVE(genelist,n=50)
#'



SIEVE<-function(genelist,n=50){
  a_mat<-data.frame(frequency=numeric(n),number.of.genes= numeric(n))
  a_mat$frequency<-c(1:n)
  a_res<-table(as.vector(genelist))
  for (i in 1:n){
    a_mat[i,2]<-length(a_res[a_res==i])
  }
  plot(a_mat,type = "h")
  geneset<-names(a_res[a_res>=n])
  return(geneset)
}


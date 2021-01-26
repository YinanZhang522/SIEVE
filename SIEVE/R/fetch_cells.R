#' Find the overlapped genes of gene sets defined by different method.
#' @name fetch_cells
#' @usage fetch_cells(expr,ratio=0.7,n=50)
#' @param expr The expression dataframe. Rows should be genes and columns should be cells.
#' @param ratio the Proportion of cell subset.The default value of ratio is set to 0.7, but can also be set to other values between 0 to 1.
#' @param n The times of subsampling.The default value of n is set to 50.

#' @return A list of random cell subsets.
#' @export

#' @examples celllist<-fetch_cells(expr,ratio=0.7,n=50)
fetch_cells<-function(expr,ratio=0.7,n=50){
  celllist<-as.data.frame(matrix(nrow=ratio*ncol(expr),ncol=n))
  for (i in 1:n){
    celllist[,i]<-sample(colnames(expr),ratio*length(colnames(expr)))
  }
  return(celllist)
}

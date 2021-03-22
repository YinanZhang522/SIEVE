\name{fetch_HVGs}
\alias{fetch_HVGs}
\title{Fetch the variable genes of cell subsets by different methods.}
\usage{
fetch_HVGs(dataframe,celllist,method=method,n=50)
}
\arguments{
\item{dataframe}{A gene expression dataframe, can be raw counts or normalized data, rows should be genes, columns should be cells.}

\item{celllist}{A matrix returned from fetch_cells}

\item{method}{the method used to fetch the variable genes,if your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap".if your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".}

\item{n}{The number of feature genes.The default value of n is set to 2000.}
}
\value{
A list of HVGs.
}
\description{
Fetch the variable genes of cell subsets by different methods.
}
\examples{
genelist<-fetch_HVGs(dataframe,celllist,method=method,n=2000)
}


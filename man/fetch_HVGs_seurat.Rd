\name{fetch_HVGs}
\alias{fetch_HVGs}
\title{Fetch the variable genes of cell subsets by different methods.}
\usage{
fetch_HVGs(object,celllist,method=method,n=50)
}
\arguments{
\item{object}{A Seurat object.}

\item{celllist}{A matrix returned from fetch_cells}

\item{method}{the method used to fetch the variable genes,can be chosen from "schs"(singleCellHaystack),"scmap","scran","ROGUE","Seurat_sct","M3Drop","Seurat_vst","Seurat_disp".}

\item{n}{The number of feature genes.The default value of n is set to 2000.}
}
\value{
A list of HVGs.
}
\description{
Fetch the variable genes of cell subsets by different methods.
}
\examples{
genelist<-fetch_HVGs(object,celllist,method=method,n=2000)
}


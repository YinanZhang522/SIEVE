\name{fetch_HVGs_seurat}
\alias{fetch_HVGs_seurat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fetch_HVGs_seurat(object, celllist, method = c("schs", "Scmap", "Scran", "ROGUE", "Seurat_sct", "M3Drop", "Seurat_vst", "Seurat_disp"), n = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{celllist}{
%%     ~~Describe \code{celllist} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, celllist, method = c("schs", "Scmap", "Scran", 
    "ROGUE", "Seurat_sct", "M3Drop", "Seurat_vst", "Seurat_disp"), 
    n = 2000) 
{
    celllist <- celllist
    X <- object
    X@meta.data$name <- rownames(X@meta.data)
    genelist <- matrix(nrow = n, ncol = ncol(celllist))
    for (i in 1:ncol(celllist)) {
        X.subset = subset(X, name \%in\% celllist[, i])
        X.subset = Seurat::NormalizeData(X.subset, verbose = FALSE)
        X.subset = Seurat::ScaleData(X.subset, verbose = FALSE)
        if (method == "ROGUE") {
            rogue_norm <- ROGUE::SE_fun(X.subset@assays$RNA@data)
            genelist[, i] = rogue_norm$Gene[c(1:n)]
        }
        else if (method == "scmap") {
            sce <- SingleCellExperiment::SingleCellExperiment(list(counts = X.subset@assays$RNA@counts))
            SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 
                1)
            rowData(sce)$feature_symbol <- rownames(sce)
            sce <- sce[!duplicated(rownames(sce)), ]
            sce <- scmap::selectFeatures(sce, suppress_plot = FALSE, 
                n_features = n)
            genelist[, i] = rowData(sce)$feature_symbol[rowData(sce)$scmap_features == 
                "TRUE"]
        }
        else if (method == "scran") {
            scran_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = X.subset@assays$RNA@counts))
            scran_sce <- scran::computeSumFactors(scran_sce)
            scran_sce <- scater::logNormCounts(scran_sce)
            stats <- scran::modelGeneVar(scran_sce)
            genelist[, i] <- scran::getTopHVGs(stats, n = n)
        }
        else if (method == "schs") {
            X.subset = Seurat::FindVariableFeatures(X.subset, 
                selection.method = "vst", nfeatures = n, verbose = FALSE)
            X.subset = Seurat::RunPCA(X.subset, npcs = 30, verbose = FALSE)
            SCHS <- singleCellHaystack::haystack(X.subset, assay = "RNA")
            genelist[, i] <- singleCellHaystack::show_result_haystack(SCHS, 
                n = n)
        }
        else if (method == "M3Drop") {
            abc_MB <- M3Drop::BrenneckeGetVariableGenes(X.subset@assays$RNA@counts, 
                fdr = 0.5)
            genelist[, i] <- abc_MB$Gene[c(1:n)]
        }
        else if (method == "Seurat_sct") {
            X.SCT <- sctransform::SCTransform(X.subset, assay = "RNA", 
                variable.features.n = n)
            genelist[, i] = X.SCT@assays$SCT@var.features
        }
        else if (method == "Seurat_vst") {
            X.subset = Seurat::FindVariableFeatures(X.subset, 
                selection.method = "vst", nfeatures = n, verbose = FALSE)
            genelist[, i] <- X.subset@assays$RNA@var.features
        }
        else if (method == "Seurat_disp") {
            X.subset = Seurat::FindVariableFeatures(X.subset, 
                selection.method = "disp", nfeatures = n, verbose = FALSE)
            genelist[, i] <- X.subset@assays$RNA@var.features
        }
    }
    return(genelist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

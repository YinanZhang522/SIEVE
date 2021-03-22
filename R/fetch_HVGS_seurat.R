#' Fetch the variable genes of cell subsets by different methods.
#' @name fetch_HVGs_seurat
#' @usage fetch_HVGs_seurat(object,celllist,method=method,n=50)
#' @param object A Seurat object.
#' @param celllist A matrix returned from fetch_cells
#' @param method the method used to fetch the variable genes,can be chosen from "schs"(singleCellHaystack),"scmap","scran","ROGUE","Seurat_sct","M3Drop","Seurat_vst","Seurat_disp".
#' @param n The number of feature genes.The default value of n is set to 2000.


#' @return A list of HVGs.
#' @export

#' @examples genelist<-fetch_HVGs_seurat(object,celllist,method=method,n=2000)
fetch_HVGs_seurat<-function(object,celllist,method=c("schs","Scmap","Scran","ROGUE","Seurat_sct","M3Drop","Seurat_vst","Seurat_disp"),n=2000){
  celllist<-celllist
  X<-object
  X@meta.data$name<-rownames(X@meta.data)
  genelist<-matrix(nrow=n,ncol=ncol(celllist))
  for (i in 1:ncol(celllist)){
    X.subset=subset(X,name%in%celllist[,i])
    X.subset=Seurat::NormalizeData(X.subset, verbose=FALSE)
    X.subset=Seurat::ScaleData(X.subset, verbose=FALSE)
    if (method=="ROGUE"){
      rogue_norm<-ROGUE::SE_fun(X.subset@assays$RNA@data)
      genelist[,i]=rogue_norm$Gene[c(1:n)]

    }else if(method=="scmap"){

      sce <- SingleCellExperiment::SingleCellExperiment(list(counts=X.subset@assays$RNA@counts))
      SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 1)
      rowData(sce)$feature_symbol <- rownames(sce)# use gene names as feature symbols
      sce <- sce[!duplicated(rownames(sce)), ]# remove features with duplicated names
      sce <- scmap::selectFeatures(sce, suppress_plot = FALSE,n_features = n)
      genelist[,i]=rowData(sce)$feature_symbol[rowData(sce)$scmap_features=="TRUE"]

    }else if(method=="scran"){
      scran_sce <- SingleCellExperiment::SingleCellExperiment(list(counts=X.subset@assays$RNA@counts))
      scran_sce <- scran::computeSumFactors(scran_sce)
      scran_sce <- scater::logNormCounts(scran_sce)
      stats <- scran::modelGeneVar(scran_sce)
      genelist[,i]<-scran::getTopHVGs(stats,n=n)

    }else if(method=="schs"){
      X.subset=Seurat::FindVariableFeatures(X.subset, selection.method="vst", nfeatures=n, verbose=FALSE)
      X.subset=Seurat::RunPCA(X.subset, npcs=30, verbose=FALSE)
      SCHS<-singleCellHaystack::haystack(X.subset,assay="RNA")
      genelist[,i]<-singleCellHaystack::show_result_haystack(SCHS, n = n)

    }else if(method=="M3Drop"){
      abc_MB<-M3Drop::BrenneckeGetVariableGenes(X.subset@assays$RNA@counts,fdr = 0.5)
      genelist[,i]<-abc_MB$Gene[c(1:n)]

    }else if(method=="Seurat_sct"){
      X.SCT<-sctransform::SCTransform(X.subset,assay="RNA",variable.features.n = n)
      genelist[,i]=X.SCT@assays$SCT@var.features

    }else if(method=="Seurat_vst"){
      X.subset=Seurat::FindVariableFeatures(X.subset, selection.method="vst", nfeatures=n, verbose=FALSE)
      genelist[,i]<-X.subset@assays$RNA@var.features

    }else if(method=="Seurat_disp"){
      X.subset=Seurat::FindVariableFeatures(X.subset, selection.method="disp", nfeatures=n, verbose=FALSE)
      genelist[,i]<-X.subset@assays$RNA@var.features
    }
  }
  return(genelist)

}

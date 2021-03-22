#' Fetch the variable genes of cell subsets by different methods.
#' @name fetch_HVGs
#' @usage fetch_HVGs(dataframe,celllist,method=method,n=50)
#' @param dataframe A gene expression dataframe, can be raw counts or normalized data, rows should be genes, columns should be cells.
#' @param celllist A data frame returned from fetch_cells
#' @param method the method used to fetch the variable genes,
#'               if your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap".
#'               if your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".
#' @param n The number of feature genes.The default value of n is set to 2000.


#' @return A list of HVGs.
#' @export

#' @examples genelist<-fatch_HVGs(data frame,celllist,method=method,n=2000)
fetch_HVGs<-function(dataframe,celllist,method=c("schs","Scmap","Scran","ROGUE","Seurat_sct","M3Drop","Seurat_vst","Seurat_disp"),n=2000){
  celllist<-celllist
  X<-dataframe
  genelist<-matrix(nrow=n,ncol=ncol(celllist))
  for (i in 1:ncol(celllist)){
    X.subset.data=X[,colnames(X)%in%celllist[,i]]

    if (method=="ROGUE"){
      rogue_norm<-ROGUE::SE_fun(X.subset.data)
      genelist[,i]=rogue_norm$Gene[c(1:n)]

    }else if(method=="scmap"){

      sce <- SingleCellExperiment::SingleCellExperiment(list(counts=X.subset.data))
      SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 1)
      rowData(sce)$feature_symbol <- rownames(sce)# use gene names as feature symbols
      sce <- sce[!duplicated(rownames(sce)), ]# remove features with duplicated names
      sce <- scmap::selectFeatures(sce, suppress_plot = FALSE,n_features = n)
      genelist[,i]=rowData(sce)$feature_symbol[rowData(sce)$scmap_features=="TRUE"]

    }else if(method=="scran"){
      scran_sce <- SingleCellExperiment::SingleCellExperiment(list(counts=X.subset.data))
      scran_sce <- scran::computeSumFactors(scran_sce)
      scran_sce <- scater::logNormCounts(scran_sce)
      stats <- scran::modelGeneVar(scran_sce)
      genelist[,i]<-scran::getTopHVGs(stats,n=n)

    }else if(method=="schs"){
      X.subset=Seurat::RunPCA(X.subset, npcs=30, verbose=FALSE)
      detection=t(X.subset.data>0)
      SCHS<-singleCellHaystack::haystack(X.subset.data,detection=detection)
      genelist[,i]<-show_result_haystack(SCHS, n = n)

    }else if(method=="M3Drop"){
      MB<-M3Drop::BrenneckeGetVariableGenes(X.subset.data,fdr = 0.5)
      genelist[,i]<-MB$Gene[c(1:n)]

    }else if(method=="Seurat_sct"){
      X.subset=Seurat::CreateSeuratObject(counts=X.subset.data, project="hspc", min.cells=50, min.features=1000,names.delim=".")
      X.subset<- Seurat::NormalizeData(X.subset)
      X.subset<- Seurat::ScaleData(X.subset)
      X.subset<-sctransform::SCTransform(X.subset,assay="RNA",variable.features.n = n)
      genelist[,i]=X.subset@assays$SCT@var.features

    }else if(method=="Seurat_vst"){
      X.subset=Seurat::CreateSeuratObject(counts=X.subset.data, project="hspc", min.cells=50, min.features=1000,names.delim=".")
      X.subset<- Seurat::NormalizeData(X.subset)
      X.subset<- Seurat::ScaleData(X.subset)
      X.subset<- Seurat::FindVariableFeatures(X.subset, selection.method = "vst", nfeatures = n)
      genes.use<-X.subset@assays$RNA@var.features

    }else if(method=="Seurat_disp"){
      X.subset=Seurat::CreateSeuratObject(counts=X.subset.data, project="hspc", min.cells=50, min.features=1000,names.delim=".")
      X.subset<- Seurat::NormalizeData(X.subset)
      X.subset<- Seurat::ScaleData(X.subset)
      X.subset<- Seurat::FindVariableFeatures(X.subset, selection.method = "disp", nfeatures = n)
      genelist[,i]<-X.subset@assays$RNA@var.features
    }
  }
  return(genelist)

}

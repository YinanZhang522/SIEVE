# Overview
   SIEVE is a strategy used for select highly Variable genes from single-cell RNA Seq data, based on random sampling for all single cells in a scRNA-seq dataset, SIEVE provides the reproducibility estimation for HVGs selection method and screens robust HVGs for the following analysis.

# Installing SIEVE
		if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
		devtools::install_github("YinanZhang522/SIEVE")

# Usage
Fetch_cells: It's a function using for fetch cell subsets with the random sampling strategy. 
Users should provide a gene expression dataframe.Rows should be genes and columns should be cells. 
		You can choose the proportion of the cell subset.To ensure the accuracy of subsequent analysis results, 
	you should consider the size and the composition of your dataset, to make sure that the proportion of the cell subset can covered every cell type, 
	the default value of ratio is set to 0.7, but can also be set to other values between 0 to 1. 
	        You can also choose the times of subsampling according to the size of dataset and the computing time. The default value of n is set to 50. 
   
		usage: fetch_cells(expr,ratio=0.7,n=50)
		param expr: The expression dataframe. Rows should be genes and columns should be cells.
		param ratio: the Proportion of cell subset.The default value of ratio is set to 0.7, but can also be set to other values between 0 to 1.
		param n: The times of subsampling.The default value of n is set to 50.  

  
Fetch_HVGs: It's a function using for fetch the variable genes of cell subsets by different methods.
		Users should provide a gene expression dataframe, both raw counts and normalized data are supported. 
		Rows should be genes and columns should be cells. 
		The feature selecting methods can be choosing from singleCellHaystack, Scmap, Scran, ROGUE, M3Drop, Seurat(vst,disp,SCTransform).
		According to the type of data, 
	if your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap".
	If your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".
				
		usage fetch_HVGs(object,celllist,method=method,n=50)
		param dataframe: A gene expression dataframe, can be raw counts or normalized data, rows should be genes, columns should be cells.
		param celllist: A matrix returned from fetch_cells
		param method: the method used to fetch the variable genes,
		              if your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap".
                              if your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".
		param n: The number of feature genes.The default value of n is set to 2000.


SIEVE: It's a function using for find the high frequency genes of gene sets defined by different method. 
		Users should provide a genelist matrix.Rows should be genes and columns should be different batches. 
		You can set the minimum number of repetitions of a gene. 

		usage: SIEVE(genelist,n=50)
		param genelist The gene sets matrix. Rows should be genes and columns should be different batches.
		param n: The minimum number of repetitions of a gene.




# Contact
	Please contact us:
	YinanZhang:zhangyinan@ihcams.ac.cn
	PengWu:wupeng1@ihcams.ac.cn

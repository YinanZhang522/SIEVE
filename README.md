# Overview
		SIEVE is a strategy used for select highly Variable genes from single-cell RNA Seq data, based on random sampling for all single cells in a scRNA-seq dataset, SIEVE provides the reproducibility estimation for HVGs selection method and screens robust HVGs for the following analysis.

# Installing SIEVE
		if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
		devtools::install_github("YinanZhang522/SIEVE")

# Usage
		Fetch_cells: Find the overlapped genes of gene sets defined by different method.
   
		usage: fetch_cells(expr,ratio=0.7,n=50)
		param expr: The expression dataframe. Rows should be genes and columns should be cells.
		param ratio: the Proportion of cell subset.The default value of ratio is set to 0.7, but can also be set to other values between 0 to 1.
		param n: The times of subsampling.The default value of n is set to 50.  

  
		Fetch_HVGs: Fetch the variable genes of cell subsets by different methods.
   
		usage fetch_HVGs(object,celllist,method=method,n=50)
		param object: A Seurat object.
		param celllist: A matrix returned from fetch_cells
		param method: the method used to fetch the variable genes,can be chosen from "schs"(singleCellHaystack),"scmap","scran","ROGUE".
		param n: The number of feature genes.The default value of n is set to 2000.


		SIEVE: Find the high frequency genes of gene sets defined by different method.

		usage: SIEVE(genelist,n=50)
		param genelist The gene sets matrix. Rows should be genes and columns should be different batches.
		param n: The minimum number of repetitions of a gene.




# Contact
Please contact us:
YinanZhang:zhangyinan@ihcams.ac.cn
PengWu:wupeng1@ihcams.ac.cn

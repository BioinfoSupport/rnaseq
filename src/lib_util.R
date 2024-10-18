
library(SummarizedExperiment)
library(tidyverse)



# Helper function to extract a data.frame of all gene-expression values
# suitable to display a heatmap
makeAssayDF <- function(x,assay_idx=1,col_meta=colnames(colData(x)),row_meta=character(0)) {
  assay(x,assay_idx) |>
    as.data.frame() |>
    rownames_to_column("row_id") |>
    pivot_longer(!row_id,names_to = "col_id") |>
    full_join(as.data.frame(colData(x)[col_meta]) |> rownames_to_column("col_id"),by="col_id") |> # Add all colData
    full_join(as.data.frame(rowData(x)[row_meta]) |> rownames_to_column("row_id"),by="row_id") # Add all rowData
}

fct_reorder_hclust <- function(.f,.x,.f2,.fun=mean) {
	h <- tibble(.f,.x,.f2) |>
		group_by(.f,.f2) |>
		summarize(.x=mean(.x)) |>
		ungroup() |>
		pivot_wider(id_cols = ".f",names_from = ".f2",values_from = ".x") |>
		column_to_rownames(".f") |>
		as.matrix() |>
		dist() |>
		hclust(method="ward.D2")
	fct_relevel(.f,h$label[h$order])
}
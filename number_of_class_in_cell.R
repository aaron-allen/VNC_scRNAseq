#!/usr/bin/env Rscript

# Aaron M. Allen, 2019.11.23

# Function to return a list of cells and the number of genes from a list of 
# genes that each cell expresses. 
# This function expects 4 imputs:
#       object = a seurat object (only tested for seurat v2)
#       gene_list = a list of gene names (for my data, this is R6.13 gene symbols)
#       accept_low =  numeric, a lower threshold for expression
#       use_raw = logical, whether to use the raw count data or the normalized data
#       use_scale = logical, whether to use the scaled data or the normalized data
#       cells_use = list, which cells to include



number_of_class_in_cell <- function(object,
                                    gene_list,
                                    accept_low = 0,
                                    cells_use = NULL,
                                    use_raw = FALSE,
                                    use_scaled = FALSE) {
    
    if (length(setdiff(gene_list, object@raw.data@Dimnames[[1]]))>0) {
        cat("The following genes were not found in the dataset:\n",
            setdiff(gene_list, object@raw.data@Dimnames[[1]]))
    }
    exp_gene_list <- intersect(gene_list, object@raw.data@Dimnames[[1]])
    exp_gene_data <- FetchData(object = object,
                               vars.all = exp_gene_list,
                               cells.use = cells_use,
                               use.raw = use_raw,
                               use.scaled = use_scaled)
    counts <- dplyr::tibble(cell_id = rownames(exp_gene_data),
                            count = rowSums(exp_gene_data>accept_low))
    return(counts)
}    

####

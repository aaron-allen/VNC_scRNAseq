#!/usr/bin/env Rscript

# Aaron M. Allen, 2019.11.23

# Function to return a list of cells and the number of genes from a list of 
# genes that each cell expresses. 
# This function expects 4 imputs:
#       object = a seurat object (only tested for seurat v2)
#       gene_list = a list of gene names (for my data, this is R6.13 gene symbols)
#       accept_low =  numeric, a lower threshold for expression
#       use_raw = logical, whether to use the raw count data or the normalized data



## Much (much, much, much, much, ...) faster version using Seurat's FetchData
## function. Now can also pass use.scaled and cells.use options.

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

######

# if (length(setdiff(
#     object@raw.data@Dimnames[[2]],
#     rownames(exp_gene_list)))>0) {
#     nonexpr_cells <- dplyr::tibble(cell_id = setdiff(object@raw.data@Dimnames[[2]],
#                                                      rownames(exp_gene_data)),
#                                    count = 0)
#     counts <- dplyr::bind_rows(counts,nonexpr_cells)
# }

######

# number_of_class_in_cell <- function(object, gene_list, accept_low, use_raw) {
#     
#     tallied <- dplyr::tibble(cell_id = object@raw.data@Dimnames[[2]],
#                              count = integer(dim(object@raw.data)[2]))
#     sub_gene_list <- intersect(gene_list, object@raw.data@Dimnames[[1]])
#     curr_progress = 0
#     for (gene in sub_gene_list) {
#         if (use_raw) {
#             temp <- WhichCells(object = object,
#                                subset.name = gene,
#                                accept.low = accept_low,
#                                use.raw = TRUE)
#         } else {
#             temp <- WhichCells(object = object,
#                                subset.name = gene,
#                                accept.low = accept_low,
#                                use.raw = FALSE)
#         }
#         if (length(temp)>0) {
#             for (ii in 1:length(temp)) {
#                 tallied$count[tallied$cell_id == temp[ii]] = 
#                     tallied$count[tallied$cell_id == temp[ii]] + 1
#             }
#         }
#         
#         # # progress bar ...
#         curr_progress = curr_progress + 1
#         curr_progress_percent = round((100*curr_progress)/length(sub_gene_list),0)
#         Sys.sleep(1)
#         print(curr_progress_percent)
#         
#         # svMisc::progress(curr_progress, progress.bar = TRUE)
#         # Sys.sleep(0.01)
#         # if (curr_progress_percent == 101) cat("Done!\n")
#         
#     }
#     return(tallied)
# }
# 
# 


###### crazy slow version - up to 8.5 times slower

# number_of_class_in_cell <- function(object, gene_list, accept_low, use_raw) {
#     tally <- integer(dim(object@raw.data)[2])
#     sub_gene_list <- intersect(gene_list, object@raw.data@Dimnames[[1]])
#     for (gene in sub_gene_list) {
#         if (use_raw) {
#             temp <- as.integer(
#                         as.logical(
#                             object@raw.data[object@raw.data@Dimnames[[1]]==gene]>accept_low))
#         } else {
#             temp <- as.integer(
#                         as.logical(
#                             object@data[object@data@Dimnames[[1]]==gene]>accept_low))
#                 }
#         tally <- tally + temp
#     }
#     tallied <- tibble(cell_id = object@raw.data@Dimnames[[2]], count = tally)
#     return(tallied)
# }

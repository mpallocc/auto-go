#' @title Read genes list
#'
#' @description Function employed for reading all the gene lists on which the enrichment will be performed.
#' @param where_results Specify the folder in which you want to save the outputs. Default is "./". Note: if you are working with RNotebook the default working directory, if not specified, is the folder in which the .Rmd file is saved.
#' @param outfolder The name assigned to the folder in which outputs are saved. Default is: "results/". NOTE: please add "/" at the end.
#' @param log2FC_threshold Threshold value for log2(Fold Change) for considering genes as differentially expressed (default = 0).
#' @param padj_threshold Threshold value for adjusted p-value significance (Defaults to 0.05).
#' @param which_list one of c("up_genes","down_genes","up_down_genes","everything"). Select a list of genes to perform the enrichment. Respectively, both up and down regulated genes (up_down_genes), only up regulated genes (up_genes), only down regulated genes (down_genes), or everything allow to load all the three kind of lists separately.
#' @param from_DE_analysis Default is TRUE, set FALSE if the lists you want to upload are not from a differential expression analysis.where_files
#' @param where_files (Default = NULL). When from_DE_analysis = T it is mandatory to provide a path to specify where the list of genes are.
#' @param files_format (Default = NULL). when from_DE_analysis = T it is mandatory to provide the extension of the list of genes to upload.
#' @export


read_gene_list <- function(gene_lists_path = "./results",
                           log2FC_threshold = 0,
                           padj_threshold = 0.05,
                           which_list = c(
                             "up_genes",
                             "down_genes",
                             "up_down_genes",
                             "everything"
                           ),
                           from_DE_analysis = TRUE,
                           files_format = ".txt") {
  if (from_DE_analysis) {
     gene_lists_files <- list.files(path = gene_lists_path, pattern = ".*genes_list_.*.txt", recursive = TRUE, full.names = TRUE)
     to_read <- gene_lists_files[grepl(pattern = paste0("thFC", log2FC_threshold, "_thPval", padj_threshold), gene_lists_files)]

  } else if (!from_DE_analysis) {
    gene_lists_files <- list.files(path = gene_lists_path, pattern = files_format, recursive = TRUE, full.names = TRUE)
    to_read <- gene_lists_files
    which_list <- "everything"
  }

  if (which_list == "up_down_genes") {
    to_read <- to_read[grepl(pattern = "/up_down_genes", to_read)]
  } else if (which_list == "up_genes") {
    to_read <- to_read[grepl(pattern = "/up_genes", to_read)]
  } else if (which_list == "down_genes") {
    to_read <- to_read[grepl(pattern = "/down_genes", to_read)]
  } else if (which_list == "everything") {
    to_read <- to_read
  }


  gene_lists <- lapply(to_read, function(x) read.table(x, header = FALSE, sep = "\n"))
  names(gene_lists) <- tools::file_path_sans_ext(to_read)

  return(gene_lists)
}

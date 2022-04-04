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

requireNamespace("autoGO-package.R")


read_gene_list <- function(where_results = "./", outfolder = "results/", log2FC_threshold=0, padj_threshold=0.05, which_list = c("up_genes","down_genes","up_down_genes","everything"), from_DE_analysis = T, where_files = NULL, files_format = NULL) {

  if (from_DE_analysis) {
    gene_lists_path <- list.files(path = paste0(where_results,outfolder), pattern = ".*genes_list_.*.txt", recursive = T)
    gene_lists_path <- paste0(where_results,outfolder,gene_lists_path)
    to_read <- gene_lists_path[grepl(pattern = paste0("thFC",log2FC_threshold,"_thPval",padj_threshold), gene_lists_path)]
  } else if (!from_DE_analysis) {
    if (is.null(where_files)) stop("Required parameter: 'where_files' for path of gene lists.")
    if (is.null(files_format)) stop("Required parameter: 'files_format' for extension of gene lists, should be like '.txt' or '.tsv', etc.")
    gene_lists_path <- list.files(path = paste0(where_files), pattern = ".txt", recursive = T)
    gene_lists_path <- paste0(where_files,gene_lists_path)
    to_read <- gene_lists_path
    which_list <- "everything"
  }

  if (which_list == "up_down_genes") {
    to_read <- to_read[grepl(pattern = "up_down_genes", to_read)]
  } else if (which_list == "up_genes") {
    to_read <- to_read[grepl(pattern = "up_genes", to_read)]
  } else if (which_list == "down_genes") {
    to_read <- to_read[grepl(pattern = "/down_genes", to_read)]
  } else if (which_list == "everything") {
    to_read <- to_read
  }

  readed <- lapply(to_read, function (x) read.table(x, header = F, sep = '\n'))
  names(readed) <- dirname(gsub(paste0(where_results,outfolder,where_files,"|.txt"), "",to_read))

  return(readed)
}

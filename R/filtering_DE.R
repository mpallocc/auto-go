#' @title Filtering DESeq2 results
#'
#' @description We could be in a position to carry out multiple filters on the results of the differential analysis, in order not to repeat all the deseq_analysis.R code which provides for the actual computation of the differential analysis, the filtering_DE.R function has been implemented to be able to filter the file(s) " * _allres.tsv " and generate all the folders and files associated with the specific filters applied.
#' @description The function automatically searches inside the folders where_results and outfolder the file(s) (See ?deseq_analysis()) "_allres.tsv" and generates folders and files in the same folders with the new filters for foldchange and pvalue respectively.
#' @param padj_threshold (Default = 0.05) Threshold value for adjusted p-value filtering.
#' @param log2FC_threshold (Default = 0) Threshold value for log2(Fold Change) filtering.
#' @param where_results (Default = "./") Folder in which the new output is written.
#' @param outfolder (Default = "./results") Name of the folder in which the new output is written.
#' @param save_excel (Default = FALSE) Write output in MS Excel file format (.xlsx).
#' @export


filtering_DE <- function(padj_threshold = 0.05,
                         log2FC_threshold = 1,
                         where_results = "./",
                         outfolder = "./results",
                         save_excel = FALSE) {
  all_res <- list.files(
    path = paste0(where_results, outfolder),
    pattern = "_allres.tsv", recursive = TRUE
  )
  all_res <- paste0(where_results, outfolder, all_res)
  readed <- lapply(all_res, function(x) read_tsv(x, col_types = cols()))
  names(readed) <- gsub(paste0(where_results, outfolder, "|\\/.*"), "", all_res)

  for (files in names(readed)) {
    data <- readed[[files]]
    filtered <- data %>%
      dplyr::filter(.data$padj < padj_threshold & abs(.data$log2FoldChange) > log2FC_threshold)

    groups_fold <- paste0(where_results, outfolder, files, "/filtered_DE", "_thFC", log2FC_threshold, "_thPval", padj_threshold)
    groups_fold_thresh_up_down <- paste0(groups_fold, "/up_down_genes")
    groups_fold_thresh_up <- paste0(groups_fold, "/up_genes")
    groups_fold_thresh_down <- paste0(groups_fold, "/down_genes")

    # saving filtered results in different folders by thresholds
    if (!dir.exists(groups_fold)) dir.create(groups_fold, recursive = T)
    write_tsv(filtered, paste0(groups_fold, "/filtered_DE_", files, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".tsv"))
    if (save_excel) openxlsx::write.xlsx(filtered, file = paste0(groups_fold, "/filtered_DE_", files, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".xlsx"), row.names = F)

    # saving gene lists
    if (!dir.exists(groups_fold_thresh_up_down)) dir.create(groups_fold_thresh_up_down, recursive = T)
    write.table(filtered$genes, paste0(groups_fold_thresh_up_down, "/up_down_genes_list_", files, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt"), quote = F, row.names = F, col.names = F)

    if (!dir.exists(groups_fold_thresh_up)) dir.create(groups_fold_thresh_up, recursive = T)
    write.table(filtered$genes[filtered$log2FoldChange > 0], paste0(groups_fold_thresh_up, "/up_genes_list_", files, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt"), quote = F, row.names = F, col.names = F)

    if (!dir.exists(groups_fold_thresh_down)) dir.create(groups_fold_thresh_down, recursive = T)
    write.table(filtered$genes[filtered$log2FoldChange < 0], paste0(groups_fold_thresh_down, "/down_genes_list_", files, "_thFC", log2FC_threshold, "_thPval", padj_threshold, ".txt"), quote = F, row.names = F, col.names = F)
  }
}

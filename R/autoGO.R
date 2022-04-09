#' @title autoGO
#'
#' @description Perform enrichment analysis on all the desired gene lists. This function take advantage of the 'enrichR' package.
#' @param list_of_genes it can be a list of dataframe (i.e. the output of read_gene_list()), a vector or a dataframe containing a list of genes or a path to the .txt file.
#' @param dbs Databases over which the enrichment will be performed, based on the enrichR libraries. Default are GO_Molecular_Function_2021, GO_Cellular_Component_2021, GO_Biological_Process_2021, KEGG_2021_Human.
#' @param my_comparison Name of the comparison the user would like to inspect.
#' @param ensembl (Default = FALSE). Set to TRUE if the provided gene list contains Ensembl IDs. A conversion to HGNC will be performed.
#' @param excel (Default = FALSE). Set to TRUE if you want to save output tables in .xlsx format.
#' @param where_results Specify the folder in which you want to save the outputs. Default is "./". Note: if you are working with RNotebook the default working directory, if not specified, is the folder in which the .Rmd file is saved.
#' @param outfolder The name to assign to the folder in which outputs are saved. Default is: "results/". NOTE: please add "/" at the end.
#' @param my_autoGO_dir where have you cloned the auto-go repository, default is your home directory "~/"
#' @export
#'
#' @import tidyverse
#' @import readr
#' @import dplyr
#' @import gdata
#' @import reshape2
#' @import circlize
#' @import DESeq2
#' @import ComplexHeatmap
#' @import enrichR
#' @import GSVA
#' @import utils
#' @import ggplot2
#' @import stringr
#' @import stats
#' @import SummarizedExperiment
#' @import textshape
#' @import ape
#' @import openxlsx
#' @import tidyr
#' @import dichromat
#' @import grid
#' @import purrr
#' @import Rcssplot
#' @import imguR
#' @import radiant.data
#' @import ggrepel
#' @import grDevices




autoGO <- function(list_of_genes, dbs = c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human"),
                   my_comparison, ensembl = F, excel = F, where_results = "./", outfolder = "results/", my_autoGO_dir = "~/") {

  if (is.data.frame(list_of_genes)) {
    list_of_genes <- list_of_genes %>% pull()
  } else if (is.character(list_of_genes) & !grepl(".txt", list_of_genes)[1]) {
    list_of_genes <- list_of_genes
  } else if (grepl(".txt", list_of_genes)) {
    list_of_genes <- read_delim(list_of_genes, delim = '\t', col_types = cols(), col_names = F) %>% pull()
  }

  if (ensembl) {
    conv_path <-  "data/conversion_ensembl_hgnc.txt"
    print(conv_path)
    all_genes_conversion <- read_delim(conv_path, delim = '\t', col_types = cols())
    list_of_genes <- as.data.frame(list_of_genes) %>% inner_join(all_genes_conversion, by=c("list_of_genes"="ensembl_gene_id")) %>%  pull()
  }

  enriched <- enrichr(list_of_genes, dbs)
  tables <- lapply(seq_along(enriched), function(i) {
    enriched[[i]] <- enriched[[i]][!is.na(enriched[[i]]$Term),]
    enriched[[i]]$`-log10(Adjusted.P.value)` <- -log10(enriched[[i]]$Adjusted.P.value)
    enriched[[i]] <- enriched[[i]][order(enriched[[i]]$Adjusted.P.value),]
  })

  if (!grepl("\\./results/", my_comparison)) {
    my_path <- paste0(where_results,outfolder,my_comparison,"/enrichment_tables/")
    if (!dir.exists(my_path)) dir.create(my_path, recursive = T)
  } else {
    my_path <- paste0(gsub("[^\\/]+$","",my_comparison),"/enrichment_tables/")
    if (!dir.exists(my_path)) dir.create(my_path, recursive=T)
  }

  invisible(lapply(seq_along(enriched), function(ind) {
    if (dim(enriched[[ind]])[1] > 0 ) {
      if(excel) xlsx::write.xlsx(enriched[[ind]], file=paste0(my_path,names(enriched)[ind], ".xlsx"),row.names = F)
      write.table(enriched[[ind]], sep="\t", quote=F, file=paste0(my_path,names(enriched)[ind], ".tsv"),row.names = F)
      }
  }))

}

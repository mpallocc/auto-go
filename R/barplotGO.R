#' @title barplotGO
#'
#' @description The function barplotGO.R implement the barplot of the first 15 enrichment terms.
#' @description For each enrichment result table a barplot is produced. Results are stored in the "enrichment_plots" subfolder for each comparison.
#' @param enrich_table Dataframe containing the enrichment results or a path to your .tsv file containing the enrichment results. Columns 'Term' and 'Adjusted.P.Value' are required.
#' @param my_comparison Name of the comparison the user would like to inspect.
#' @param where_results Specify the folder in which you want to save outputs. (Default = "./"). Note: if you are working with R Notebooks the default working directory (if not specified) is the folder in which the .Rmd is saved.
#' @param outfolder The name to assign to the folder for output saving. (Default = "results/"). NOTE: please add "/" at the end.
#' @export


barplotGO <- function(enrich_table,  my_comparison = NULL, where_results = "./", outfolder = "results/") {

  if (is.data.frame(enrich_table)) {
    enrich_table <- enrich_table
  } else if (grepl(".tsv", enrich_table)) {
    enrich_table_path <- enrich_table
    enrich_table <- read_delim(enrich_table, delim = "\t", col_types = cols())
  }


  if (is.null(my_comparison)) {
      my_analysis <- str_match(enrich_table_path, pattern = paste0(where_results,outfolder,"(.*?)\\/"))[2]
      dbs <- gsub(".*\\/|\\.tsv", "", enrich_table_path)
      path_to_save <- paste0(gsub("_tables.*","", enrich_table_path),"_plots")
  } else if (!is.null(my_comparison)) {
    my_analysis <- str_match(my_comparison, pattern = paste0(where_results,outfolder,"(.*?)\\/"))[2]
    if (grepl("\\/", my_comparison)) {
      dbs <- gsub(".*\\/", "", my_comparison)
      path_to_save <- paste0(gsub("_tables.*","", my_comparison),"_plots/")
    } else if (!grepl("\\/", my_comparison)){
      dbs <- NULL
      path_to_save <- paste0("./",my_analysis,"/enrichment_plots/")
    }
  }

  if (grepl("down_genes", path_to_save)) {
    title <- paste0(gsub("_"," ",dbs), " for Down Regulated Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else if (grepl("up_genes", path_to_save)) {
    title <- paste0(gsub("_"," ",dbs), " for Up Regulated Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else if (grepl("up_down_genes", path_to_save)) {
    title <- paste0(gsub("_"," ",dbs), " for all DE Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else {
    title <- ifelse(is.na(gsub("_", " ", dbs)), "", gsub("_", " ", dbs))
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  }

  enrich_table <- enrich_table %>%
    arrange(Adjusted.P.value, Term) %>%
    dplyr::slice(1:15) %>%
    mutate(Term = gsub("\\(GO.*","",Term),
           `-log10(Adjusted.P.value)` = -log10(Adjusted.P.value))

  ggplot(data=enrich_table, aes(x=reorder(.data$Term, `-log10(Adjusted.P.value)`), y=`-log10(Adjusted.P.value)`)) +
    geom_bar(stat="identity",fill="#91bbdb", width = 0.5) +
    theme_minimal() +
    ggtitle(title, subtitle = subtitle) +
    geom_hline(yintercept = -log10(0.05), size = 1, colour = "#e09696", linetype="longdash") +
    scale_y_continuous(expand = c(0,0))+
    coord_flip() +
    labs(y = "-log10(adjp_val)", x = "") +
    theme(title = element_text(size = 23, ), plot.background = element_rect(fill = "#ffffff"),
          axis.text=element_text(size=18), axis.title = element_text(size = 20),
          axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))

  if (!dir.exists(path_to_save)) dir.create(path_to_save, recursive=T)

  ggsave(filename=paste0(path_to_save,"barplotGO_",dbs,".png"), plot=last_plot(), width = unit(25,'cm'), height = unit(10,'cm'))

}

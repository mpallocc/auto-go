#' @title lolliGO
#'
#' @description The function lolliGO.R implement the lollipop plot of the first 20 enrichment terms.
#' @description For each enrichment result table a lollipop plot is produced. Results are stored in the "enrichment_plots" subfolder for each comparison.
#' @param enrich_table Dataframe containing the enrichment results or a path to your .tsv file containing the enrichment results. Columns 'Term' and 'Adjusted.P.Value' are required. It has to contain also the column Overlap which indicates the number of genes in our list enriched over the total genes of that term (format n/N).
#' @param my_comparison Name of the comparison the user would like to inspect.
#' @param where_results Specify the folder in which you want to save outputs. (Default = "./"). Note: if you are working with R Notebooks the default working directory (if not specified) is the folder in which the .Rmd is saved.
#' @param outfolder The name to assign to the folder for output saving. (Default = "results/"). NOTE: please add "/" at the end.
#' @export


lolliGO <- function(enrich_table,
                    my_comparison = NULL,
                    where_results = "./",
                    outfolder = "results/") {
  if (is.data.frame(enrich_table)) {
    enrich_table <- enrich_table
  } else if (grepl(".tsv", enrich_table)) {
    enrich_table_path <- enrich_table
    enrich_table <- read_delim(enrich_table, delim = "\t", col_types = cols())
  }

  pattern <- paste0(where_results, outfolder, "(.*?)\\/")
  if (is.null(my_comparison)) {
    my_analysis <- str_match(enrich_table_path,
      pattern = pattern
    )[2]
    dbs <- gsub(".*\\/|\\.tsv", "", enrich_table_path)
    path_to_save <- paste0(gsub("_tables.*", "", enrich_table_path), "_plots")
  } else if (!is.null(my_comparison)) {
    my_analysis <- str_match(my_comparison, pattern = pattern)[2]
    if (grepl("\\/", my_comparison)) {
      dbs <- gsub(".*\\/", "", my_comparison)
      path_to_save <- paste0(gsub("_tables.*", "", my_comparison), "_plots/")
    } else if (!grepl("\\/", my_comparison)) {
      dbs <- NULL
      path_to_save <- paste0("./", my_analysis, "/enrichment_plots/")
    }
  }

  if (grepl("/down_genes", path_to_save)) {
    title <- paste0(gsub("_", " ", dbs), " for Down Regulated Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else if (grepl("up_genes", path_to_save)) {
    title <- paste0(gsub("_", " ", dbs), " for Up Regulated Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else if (grepl("up_down_genes", path_to_save)) {
    title <- paste0(gsub("_", " ", dbs), " for all DE Genes")
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  } else {
    title <- ifelse(is.na(gsub("_", " ", dbs)), "", gsub("_", " ", dbs))
    subtitle <- ifelse(is.na(gsub("_", " ", my_analysis)), "", gsub("_", " ", my_analysis))
  }

  enrich_table <- enrich_table %>%
    dplyr::arrange(.data$Adjusted.P.value, .data$Term) %>%
    dplyr::slice(1:20) %>%
    tidyr::extract(.data$Overlap, into = c("gene_counts", "gene_total"), regex = "([0-9]+)\\/([0-9]+)") %>%
    type_convert(col_types = cols(gene_counts = col_double(), gene_total = col_double())) %>%
    dplyr::mutate(
      Term = gsub("\\(GO.*", "", .data$Term),
      `-log10(Adjusted.P.value)` = -log10(.data$`Adjusted.P.value`),
      percent = round(.data$gene_counts / .data$gene_total, digits = 2),
      percent = ifelse(.data$percent > 0.4, 0.4, .data$percent)
    )

  breaks <- round(seq(min(enrich_table$gene_counts), max(enrich_table$gene_counts), length.out = 6))

  ggplot(enrich_table, aes(x = .data$`-log10(Adjusted.P.value)`, reorder(.data$Term, .data$`Adjusted.P.value`))) +
    ggtitle(label = title, subtitle = subtitle) +
    geom_segment(aes(xend = 0, yend = .data$Term)) +
    geom_point(aes(color = .data$percent, size = .data$gene_counts)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE), option = "plasma", breaks = seq(0, 0.4, 0.1), limits = c(0, 0.4), labels = c("0 %", "10 %", "20 %", "30 %", "> 40 %")) +
    scale_size_continuous(range = c(5, 12), breaks = breaks) +
    theme_minimal() +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = "-log10(adjp_val)", y = "") +
    geom_vline(xintercept = -log10(0.05), size = 1, colour = "#e09696", linetype = "longdash") +
    theme(
      title = element_text(size = 23), plot.background = element_rect(fill = "#ffffff"),
      axis.text = element_text(size = 18), axis.title = element_text(size = 20),
      axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0))
    ) +
    guides(
      colour = guide_colourbar(title = "Percentage", order = 1, title.position = "top", title.theme = element_text(size = 15), label.vjust = 0.5, label.theme = element_text(size = 12), ticks.colour = "black"),
      size = guide_legend(title = "Counts", order = 2, title.position = "top", title.theme = element_text(size = 15), reverse = T, label.theme = element_text(size = 12))
    )

  if (!dir.exists(path_to_save)) dir.create(path_to_save, recursive = T)
  ggsave(filename = paste0(path_to_save, "lolliGO_", dbs, ".png"), plot = last_plot(), width = unit(20, "cm"), height = unit(10, "cm"))
}

#' @title HeatmapGO
#'
#' @description If the analysis has been performed on more conditions it is interest to have a look at the difference in the enrichment results between the groups. This can be performed by heatmapGO().
#' @description The function automatically reads all the enrichment results of the chosen database. A heatmap is produced for each database, all the terms are merged together and a filter is applied as follows: only terms with a significant pvalue (i.e. less than padj_threshold) in at least one comparison will be retained and plotted. These plots will be saved in the "Comparison_Heatmap" folder. In order to have readable plots, if many terms are enriched for a database several images will be created (indexed _1, _2, ...).
#' @param lib Database of choice to plot the heatmap. It has to be one for which the enrichment analysis has been performed.
#' @param where_results Specify the folder in which you want to save outputs. (Default = "./"). Note: if you are working with R Notebooks the default working directory (if not specified) is the folder in which the .Rmd is saved.
#' @param outfolder The name to assign to the folder for output saving. (Default = "./results"). NOTE: please add "/" at the end.
#' @param log2FC_threshold Threshold value for log2(Fold Change) for considering genes as differentially expressed (Default = 0).
#' @param padj_threshold Threshold value for adjusted p-value significance (Defaults to 0.05).
#' @param which_list One of c("up_genes", "down_genes","up_down_genes", "not_from_DE"): select data to plot. Respectively, only up regulated genes (up_genes), only down regulated genes ("down_genes"), enrichment on both up and down regulated genes (up_down_genes) or select "not_from_DE" if the enrichment will be made on a list of genes that does not come from a differential expression analysis.
#' @export

# TODO: this function should be told:
# * where enrichment results are
# ** we can provide no `lib` (it will produce heatmaps from all enrichment DBs)
# ** user can provide a list of `libs` and in that case we produce all the heatmaps?
# * where to save the heatmaps: `where_results` & `outfolder`

heatmapGO <- function(lib,
                      where_results = "./",
                      outfolder = "./results",
                      log2FC_threshold = 0,
                      padj_threshold = 0.05,
                      which_list = c(
                        "up_genes",
                        "down_genes",
                        "up_down_genes",
                        "not_from_DE"
                      )) {
  x <- list.files(
    pattern = paste0(lib, ".tsv"),
    path = paste0(where_results, outfolder),
    recursive = TRUE,
    all.files = TRUE
  )
  x <- paste0(where_results, outfolder, x)

  if (which_list != "not_from_DE") {
    # FIXME: change column names to the new, correct ones
    to_read <- x[grepl(pattern = paste0("thFC", log2FC_threshold, "_thPval", padj_threshold), x)]
  }

  if (which_list == "up_down_genes") {
    to_read <- to_read[grepl(pattern = "/up_down_genes", to_read)]
    title <- " for all DE Genes"
  } else if (which_list == "up_genes") {
    to_read <- to_read[grepl(pattern = "/up_genes", to_read)]
    title <- " for Up Regulated Genes"
  } else if (which_list == "down_genes") {
    to_read <- to_read[grepl(pattern = "/down_genes", to_read)]
    title <- " for Down Regulated Genes"
  } else if (which_list == "not_from_DE") {
    to_read <- x
    title <- ""
  }

  dd <- lapply(to_read, function(x) {
    split_path <- unlist(strsplit(x, "\\/"))
    vs_path <- split_path[grep("_vs_", split_path)]
    my_comp <- vs_path[which(min(nchar(vs_path[grep("_vs_", vs_path)])) == nchar(vs_path[grep("_vs_", vs_path)]))]
    if (identical(my_comp, character(0))) my_comp <- str_match(string = x, pattern = paste0(outfolder, "(.*)/enrichment"))[2]
    out <- read_delim(x, "\t", col_types = cols()) %>%
      dplyr::mutate(
        Comparison = my_comp,
        Comparison = ifelse(is.na(.data$Comparison), str_match(x, pattern = ("results.*\\/(.*)\\/enrichment"))[[2]], .data$Comparison)
      ) %>%
      dplyr::select(.data$Term, .data$Adjusted.P.value, .data$Comparison) %>%
      tidyr::pivot_wider(names_from = .data$Comparison, values_from = .data$Adjusted.P.value)
    return(out)
  })

  complete_table <- purrr::reduce(dd, dplyr::full_join, by = "Term") %>%
    dplyr::mutate(dplyr::across(is.na, tidyr::replace_na, 1)) %>%
    dplyr::rename_with(~ gsub("up_genes/|down_genes/", "", .x)) %>%
    dplyr::filter_all(dplyr::any_vars(.data <= padj_threshold)) %>%
    dplyr::mutate(Term = gsub("\\(GO.*", "", .data$Term)) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ (-1 * log10(.x)))) %>%
    textshape::column_to_rownames(loc = "Term")


  pretty_labels <- function(data) {
    labels <- gsub(" $", "", rownames(data))
    for (row in seq_along(labels)) {
      a <- which(strsplit(labels[row], "")[[1]] == " ")
      if (length(a) > 6) {
        substr(labels[row], a[7], a[7]) <- "\n"
      }
      if (!grepl("\n", labels[row]) & nchar(labels[row]) > 60) {
        substr(labels[row], a[length(a) - 1], a[length(a) - 1]) <- "\n"
      }
    }
    return(labels)
  }

  rownames(complete_table) <- pretty_labels(complete_table)

  plots <- function(data) {
    data <- as.matrix(data)

    col_fun <- RColorBrewer::brewer.pal(9, "YlGnBu")

    plot_name <- paste0("Multi-DE Ontology Heatmap for ", gsub("_", " ", lib), title)

    Var <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(ncol(data)), colnames(data))

    anno <- HeatmapAnnotation(
      Comparisons = colnames(data), gp = gpar(col = "black"), col = list(Comparisons = Var),
      annotation_legend_param = list(nrow = 3)
    )

    if (ncol(data) > 5) {
      width_col <- ncol(data)
    } else {
      width_col <- 6
    }

    h1 <- suppressWarnings(draw(Heatmap(data,
      name = "-log10(Adj. P-value + 1)", col = col_fun,
      row_names_gp = gpar(fontsize = 10),
      show_row_dend = F, show_column_names = F,
      show_column_dend = F, row_gap = unit(55, "cm"),
      heatmap_width = unit(2, "cm") * width_col,
      bottom_annotation = anno,
      heatmap_legend_param = list(
        legend_direction = "vertical",
        title_position = "leftcenter-rot",
        legend_height = unit(3, "cm")
      ),
      column_title = plot_name[1], rect_gp = gpar(col = "white", lwd = 1),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (data[i, j] > -log10(padj_threshold)) grid.text(sprintf("%.2f", data[i, j]), x, y, gp = gpar(fontsize = 10, col = "#e35b5b", fontface = "bold"))
      }
    ),
    heatmap_legend_side = "left", annotation_legend_side = "bottom"
    ))
  }

  if (which_list != "not_from_DE") path_save <- paste0(where_results, outfolder, "ComparisonHeatmap/", which_list)
  if (which_list == "not_from_DE") path_save <- paste0(where_results, outfolder, "ComparisonHeatmap/")

  if (!dir.exists(path_save)) dir.create(path_save, recursive = T)

  name_save <- paste0("/", lib, ".png")

  saving <- function(data, path_save) {
    if (nrow(data) > 35) {
      n_df <- round(nrow(data) / 35)
      data$group <- 1:nrow(data) %% n_df + 1
      data_list <- split(data, data$group)
      data_list <- purrr::map(data_list, ~ (.x %>% dplyr::select(-group)))
      for (ind in seq_along(data_list)) {
        png(paste0(path_save, gsub(".png", "", name_save), "_", ind, ".png"), width = 4000, height = 3500, res = 300)
        plots(data_list[[ind]])
        dev.off()
      }
    } else if (nrow(data) < 15) {
      png(paste0(path_save, name_save), width = 4000, height = 2000, res = 300)
      plots(data)
      dev.off()
    } else if (nrow(data) < 20) {
      png(paste0(path_save, name_save), width = 4000, height = 2500, res = 300)
      plots(data)
      dev.off()
    } else if (nrow(data) < 25) {
      png(paste0(path_save, name_save), width = 4000, height = 3000, res = 300)
      plots(data)
      dev.off()
    } else if (nrow(data) <= 35) {
      png(paste0(path_save, name_save), width = 4000, height = 3500, res = 300)
      plots(data)
      dev.off()
    }
  }

  saving(complete_table, path_save)
}

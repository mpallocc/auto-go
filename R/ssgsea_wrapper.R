#' @title ssGSEA
#'
#' @description Single-sample Gene Set Enrichment Analysis. This kind of analysis is recommended when there are too few samples. The function will generate an enrichment score for each sample and associated visualizations.
#' @description The function ssGSEA is implemented on top of GSVA package. It will produce an heatmap and a violin plot. If chosen, it will produce also the enrichment score table. All the necessary data are already provided with the autoGO package (as the MSigDB gene sets). It is possible to choose whether to perform ssGSEA with all the MSigDB (database) or with a sub-group.
#' @param norm_data Path of the normalized matrix. By default it is the "deseq_vst_data.txt" computed with the function "deseq_analysis()", but other normalized counts matrix (like TPMs) can be passed to the function. Requirements: first column gene names.
#' @param MSigDB_names One of c("hgnc","entrez"). The notation for gene names in the matrix.
#' @param which_gene_set Default = NULL, all the gene sets are considered. The user can specify one or more gene_sets (c1,c2,c3, c4, c5, c6, c7, c8, h ).
#' @param write_enrich_tables Default = FALSE. Set to TRUE to save the matrix with enrichment scores.
#' @param group A matrix for the annotation on the heatmap and for grouping in the statistical analysis. Can be the path to an annotation matrix or dataframe in your global environment.
#' @param my_autoGO_dir If cloned from gitlab the path to your auto-go repository, default is home directory "~/".
#' @param where_results Specify the folder in which you want to save outputs. (Default = "./"). Note: if you are working with R Notebooks the default working directory (if not specified) is the folder in which the .Rmd is saved.
#' @param outfolder The name to assign to the folder for output saving. (Default = "ssgsea/"). NOTE: please add "/" at the end.
#' @param full_names Default = FALSE, Terms full names are codified for visualization purposes. Set to TRUE to plot full names instead.
#' @param tpm_norm (Default = FALSE). Set to TRUE to perform a TPM normalization on counts matrix before the analysis.
#' @param ensembl (Default = FALSE). Set to TRUE to convert gene names from ENSEMBL to HGNC.
#' @export


ssgsea_wrapper <- function(norm_data = "results/deseq_vst_data.txt", MSigDB_names = c("hgnc","entrez"), which_gene_set = NULL, write_enrich_tables = F, group = NULL, my_autoGO_dir = "~/",  where_results = "./", outfolder = "ssgsea/", full_names = F, tpm_norm = F, ensembl = F) {

  if (grepl(".tsv", norm_data)[1]) {
    norm_data <- read_tsv(norm_data, col_types = cols())
  } else if (grepl(".csv", norm_data)[1]) {
    norm_data <- read_csv(norm_data, col_types = cols())
  } else if (grepl(".rds", norm_data)[1]) {
    norm_data <- readRDS(norm_data)
  } else if (grepl(".txt", norm_data)[1]) {
    norm_data <- read_delim(norm_data, delim = "\t", col_types = cols())
  } else if (is.data.frame(norm_data)) {
    norm_data <- norm_data
  } else {
    stop("Provide a file .tsv, .csv, .rds or .txt tab-separated")
  }

  norm_data <- as.data.frame(norm_data)
  rownames(norm_data) <- norm_data[,1]
  norm_data[,1] <- NULL

  if (ensembl) {
    conversion <- conversion_ensembl
    all_genes_conversion <- conversion %>%
      textshape::column_to_rownames(loc = "ensembl_gene_id")

    norm_data <- merge(norm_data, all_genes_conversion, by = 0) %>%
      select(-Row.names)

    if(sum(duplicated(norm_data$external_gene_name)) != 0) warning("Conversion from ENSEMBL to HGNC has produced ", sum(duplicated(norm_data$external_gene_name)), " duplicated gene names. They are going to be filtered.")

    norm_data <- norm_data %>%
      filter(!duplicated(external_gene_name)) %>%
      column_to_rownames(var = "external_gene_name")
  }

  if (tpm_norm) {
    #gene_length <-  read.table(paste0(my_autoGO_dir,"auto-go/data/gene_length.txt"), sep = "\t", header = T)

    gene_length <- gene_length %>%
      filter(external_gene_name %in% rownames(norm_data))
    norm_data <- norm_data %>%
      filter(rownames(norm_data) %in% gene_length$external_gene_name)

    if (!all(gene_length$external_gene_name == rownames(norm_data))) gene_length <- gene_length[match(rownames(norm_data), gene_length$external_gene_name),]
    gene_length <- gene_length %>% select(width) %>% pull()

    tpm <- function(counts,len) {
      x <- counts/len
      return(t(t(x)*1e6/colSums(x)))
    }
    norm_data <- tpm(norm_data,gene_length) %>% as.data.frame()
  }

  #read MSignDB genesets
  custom_file <- function(string) {
    string <- string %>% strsplit("\t")
    names(string) <- sapply(string, function(x) x[1])
    string <-lapply(string, function(x) x[-c(1,2)])
    string <- stack(string)
    return(string)
  }

  if (MSigDB_names == "hgnc") {
    gene_sets <- list.files(path = paste0(my_autoGO_dir,"auto-go/data/MSigDB"), pattern = "symbols.gmt", recursive = T, all.files = T)
    gene_sets <- paste0(my_autoGO_dir,"auto-go/data/MSigDB/",gene_sets)
    sets <- lapply(gene_sets, readLines)
    names(sets) <- gsub(".*MSigDB/| |\\.all.*","",gene_sets)
    sets <- lapply(names(sets),function(x) custom_file(sets[[x]]))
    names(sets) <- gsub(".*MSigDB/| |\\.all.*","",gene_sets)
  } else if (MSigDB_names == "entrez") {
    gene_sets <- list.files(path = paste0(my_autoGO_dir,"/auto-go/data/MSigDB"), pattern = "entrez.gmt", recursive = T, all.files = T)
    gene_sets <- paste0(my_autoGO_dir,"auto-go/data/MSigDB/",gene_sets)
    sets <- lapply(gene_sets, readLines)
    names(sets) <- gsub(".*MSigDB/| |\\.all.*","",gene_sets)
    sets <- lapply(names(sets),function(x) custom_file(sets[[x]]))
    names(sets) <- gsub(".*MSigDB/| |\\.all.*","",gene_sets)
  } else {
    stop("MSigDB_names parameter should be 'hgnc' or 'entrez'.")
  }

  if(!dir.exists(paste0(where_results,outfolder))) dir.create(paste0(where_results,outfolder), recursive=T)

  if(!is.null(which_gene_set)) {
    if (sum(which_gene_set %in% c("c1","c2","c3","c4","c5","c6","c7","c8","h")) != length(which_gene_set)) {
      stop("which_gene_set must be one (ore more) of 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'h'.")
    }
    sets <- sets[which_gene_set]
  }

  for (gs in names(sets)) {
    assign(gs, split(sets[[gs]], sets[[gs]]$ind), envir = .GlobalEnv)
    assign(gs, lapply(get(gs), function(x) x<-x[,1]), envir = .GlobalEnv)
    suppressWarnings(assign(paste0("ssgsea_", gs), gsva(as.matrix(norm_data), get(gs), verbose = F, method = "ssgsea"), envir = .GlobalEnv))
  }

  variables <- names(as.list(.GlobalEnv))
  variables <- variables[grepl("ssgsea",variables) & !grepl("wrapper", variables) & grepl(paste(which_gene_set,collapse = "|"), variables)]

  if (write_enrich_tables) {
    if(!dir.exists(paste0(where_results,outfolder,"tables/"))) dir.create(paste0(where_results,outfolder,"tables/"), recursive=T)
    for (var in variables) {
      write.table(get(var), paste0(where_results,outfolder,"tables/",var, "_EnrichmentScore.tsv"), sep = "\t", quote = F, col.names = NA, row.names = T)
    }
  }

  if (!is.null(group)){
    if (!is.data.frame(group)[1] & grepl(".txt", group)[1]) {
      group <- read_delim(group, delim = '\t', col_types = cols())
    } else if (is.data.frame(group)) {
      group <- group
      } else {
      stop("Provide a data.frame or a .txt file for group variable.")
    }

    group_class <- split.data.frame(group, group[[2]])

    for (gs in variables) {
      data <- get(gs)
      group1 <- data[,group_class[[1]][[1]]]
      group2 <- data[,group_class[[2]][[1]]]

      table <- data.frame()
      for (x in rownames(data)) {
        w_res <- wilcox.test(as.numeric(group1[x,]),as.numeric(group2[x,]))
        t_res <- t.test(as.numeric(group1[x,]),as.numeric(group2[x,]))
        tb <- data.frame(GeneSet=x, wilcoxon=w_res$p.value, t_test=t_res$p.value,
                         mean_group1 = t_res$estimate[["mean of x"]],
                         mean_group2 = t_res$estimate[["mean of y"]])
        table <- rbind.data.frame(table, tb)
      }
      table$w_adj <- p.adjust(table$wilcoxon, method = "BH")
      table$t_adj <- p.adjust(table$t_test, method = "BH")

      colnames(table)[grep("group1", colnames(table))] <- paste0("mean_",names(group_class)[1])
      colnames(table)[grep("group2", colnames(table))] <- paste0("mean_",names(group_class)[2])

      table <- table %>% arrange(w_adj)

      write.table(table, paste0(where_results,outfolder,"tables/",gsub("ssgsea_","",gs),"_stats_",
                                names(group_class)[1], "_vs_", names(group_class)[2],".txt"),
                  quote = F, sep = "\t", row.names = F)
      sign <- table$GeneSet[(table$wilcoxon <= 0.05 | table$t_test <= 0.05)]
      if (length(sign) >= 20) {
        #sign <- table$GeneSet[(table$wilcoxon <= 0.01 | table$t_test <= 0.01)]
        sign <- table$GeneSet[1:20]
      }

      group1 <- group1 %>%
        as.data.frame() %>%
        rownames_to_column(var = "geneset") %>%
        pivot_longer(cols = -geneset,names_to = "pat") %>%
        mutate(group = names(group_class)[1])

      group2 <- group2 %>%
        as.data.frame() %>%
        rownames_to_column(var = "geneset") %>%
        pivot_longer(cols = -geneset,names_to = "pat") %>%
        mutate(group = names(group_class)[2])

      df <- rbind.data.frame(group1,group2) %>%
        filter(geneset %in% sign)

      if(!dir.exists(paste0(where_results,outfolder,"plots/"))) dir.create(paste0(where_results,outfolder,"plots/"), recursive=T)

      if (!full_names) {
        cod <- data.frame(Term = unique(df$geneset), Term_ID = paste0("TERM_",1:length(unique(df$geneset))))
        write.table(cod,paste0(where_results,outfolder,"plots/codified_term_",gsub("ssgsea_","",gs), ".txt"), sep = "\t", quote = F, row.names = F)
        df <- df %>%
          left_join(cod, by = c("geneset"="Term")) %>%
          mutate(Term_ID = factor(Term_ID ,levels = paste0("TERM_",1:length(unique(df$geneset)))))
      } else if (full_names) {
        df <- df %>%
          mutate(Term_ID = geneset) %>%
          select(-geneset)
      }

      mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(df$Term_ID)))
      p<-ggplot(data = df, aes(x = Term_ID, y = value)) +
        geom_violin(aes(fill=Term_ID), show.legend = F, trim = F, scale = "width")+
        stat_summary(fun=median, geom="point", size=1, color="black", shape=18) +
        theme_bw()+
        labs(x="",y="Enrichment Score",title = paste0("Distribution of significative genesets for ", toupper(gsub("ssgsea_","",gs))))+
        theme(legend.position = "top", legend.margin = margin(0,0,0,0,"lines"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = mycolors)
      png(paste0(where_results,outfolder,"plots/","distrib_",gsub("ssgsea_","",gs),".png"), width = 4000, height = 2500, res = 300)
      print(p)
      dev.off()

      #heatmap
      group <- group %>% arrange(1)

      if (!full_names) {
        filt <- data[sign,group[[1]]] %>%
          as.data.frame() %>%
          rownames_to_column(var = "geneset") %>%
          left_join(cod, by = c("geneset"="Term")) %>%
          mutate(Term_ID = factor(Term_ID ,levels = paste0("TERM_",1:length(unique(df$geneset))))) %>%
          column_to_rownames(var = "Term_ID") %>%
          select(-geneset)
      } else if (full_names) {
        filt <- data[sign,group[[1]]] %>%
          as.data.frame()
      }

      z_filt <- t(scale(t(filt)))
      head(z_filt)

      Var <- setNames(c("#446455","#FDD262"), unique(group[[2]]))

      ha <- HeatmapAnnotation(df = group[[2]], show_annotation_name = F,
                              annotation_legend_param = list(df = list(title = " ")),
                              col = list(df = Var))
      ht<-Heatmap(z_filt, top_annotation = ha,
                       column_title = paste0("Distribution of significative genesets for ",gsub("ssgsea_","",gs)),
                       col = RColorBrewer::brewer.pal(9, "PuRd"),
                       column_names_rot = 45, heatmap_legend_param = list(title_position='leftcenter-rot'),
                       name = "zscore(ES)", show_row_dend = F,rect_gp = gpar(col = "white", lwd = 0.5))

      png(paste0(where_results,outfolder,"plots/","heatmap_",gsub("ssgsea_","",gs),".png"), width = 3000, height = 3000, res = 300)
      draw(ht, heatmap_legend_side='left', annotation_legend_side = "left", merge_legends = T)
      dev.off()
    }
  } else {
    print("Plots will be generated with the top 15 TERMS ordered by mean for sample.")
    for (gs in variables){
      data <- get(gs)
      ordered <- apply(data, 1, mean) %>% as.data.frame()
      colnames(ordered) <- "Mean"
      ordered <- ordered %>% arrange(desc(Mean))

      data <- data %>%
        as.data.frame() %>%
        rownames_to_column(var = "geneset") %>%
        filter(geneset %in% rownames(ordered)[1:15])

      if(!dir.exists(paste0(where_results,outfolder,"plots/"))) dir.create(paste0(where_results,outfolder,"plots/"), recursive=T)

      if (!full_names) {
        cod <- data.frame(Term = unique(data$geneset), Term_ID = paste0("TERM_",1:length(unique(data$geneset))))
        write.table(cod,paste0(where_results,outfolder,"plots/codified_term_",gsub("ssgsea_","",gs), ".txt"), sep = "\t", quote = F, row.names = F)
        df <- data %>%
          left_join(cod, by = c("geneset"="Term")) %>%
          mutate(Term_ID = factor(Term_ID ,levels = paste0("TERM_",1:length(unique(data$geneset))))) %>%
          select(-geneset) %>%
          pivot_longer(cols = -Term_ID)
      } else if (full_names) {
        df <- data %>%
          mutate(Term_ID = geneset) %>%
          select(-geneset) %>%
          pivot_longer(cols = -Term_ID)
      }

      mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(df$Term_ID)))
      p<-ggplot(data = df, aes(x = Term_ID, y = value)) +
        geom_violin(aes(fill=Term_ID), show.legend = F, trim = F, scale = "width")+
        stat_summary(fun=median, geom="point", size=1, color="black", shape=18) +
        theme_bw()+
        labs(x="",y="Enrichment Score",title = paste0("Distribution of significative genesets for ", toupper(gsub("ssgsea_","",gs))))+
        theme(legend.position = "top", legend.margin = margin(0,0,0,0,"lines"),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = mycolors)
      png(paste0(where_results, outfolder,"plots/distrib_",gsub("ssgsea_","",gs),".png"), width = 4000, height = 2500, res = 300)
      print(p)
      dev.off()



      if (!full_names) {
        data <- data %>%
          as.data.frame() %>%
          left_join(cod, by = c("geneset"="Term")) %>%
          mutate(Term_ID = factor(Term_ID ,levels = paste0("TERM_",1:length(unique(data$geneset))))) %>%
          column_to_rownames(var = "Term_ID") %>%
          select(-geneset)
      } else if (full_names) {
        data <- data %>%
          as.data.frame() %>%
          column_to_rownames(var = "geneset")
      }

      z_data <- t(scale(t(data)))


      if (grepl("HALLMARK_", rownames(z_data)[1])) {
        rownames(z_data) <- gsub("HALLMARK_", "", rownames(z_data), ignore.case = T)
      }


      ht<-Heatmap(z_data,
                       column_title = paste0("Distribution of significative genesets for ",gsub("ssgsea_","",gs)),
                       col = RColorBrewer::brewer.pal(9, "PuRd"),
                       column_names_rot = 45, heatmap_legend_param = list(title_position='leftcenter-rot'),
                       name = "zscore(ES)", show_row_dend = F,rect_gp = gpar(col = "white", lwd = 0.5))

      png(paste0(where_results,outfolder,"plots/heatmap_",gsub("ssgsea_","",gs),".png"), width = 3000, height = 3000, res = 300)
      draw(ht, heatmap_legend_side='left')
      dev.off()
    }
  }
}


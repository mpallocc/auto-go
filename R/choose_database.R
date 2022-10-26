#' @title choose_database
#'
#' @description It allows the user to choose the databases on which to perform the enrichment analysis. Either it returns all the possible databases or a subset of them.
#' @param dbs_search (Default = NULL), it is a string pattern based on the name of certain databases in order to return only some databases of interest
#' @export



choose_database <- function (dbs_search = NULL) {

  #setEnrichrSite("Enrichr")
  dbs_table <- listEnrichrDbs()
  dbs <- dbs_table$libraryName
  #return(dbs)

  if (!is.null(dbs_search)) {

    index <- grepl(dbs_search, dbs)
    correlated_dbs_list <- dbs[index]
    return(correlated_dbs_list)

  }

  else {
    return (dbs)
  }

}

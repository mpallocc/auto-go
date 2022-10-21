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
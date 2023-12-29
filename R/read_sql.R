#' Read SQL script
#'
#' @description
#' Converts .sql script into a text string
#'
#' @param filepath file path of .sql script
#'
#' @export
#'

read_sql <- function(filepath) {
  con = file(filepath, "r")
  sql.string <- ""

  while (TRUE){
    line <- readLines(con, n = 1)

    if ( length(line) == 0 ){
      break
    }

    line <- gsub("\\t", " ", line)

    if(grepl("--",line) == TRUE){
      line <- paste(sub("--","/*",line),"*/")
    }

    sql.string <- paste(sql.string, line)
  }

  close(con)
  return(sql.string)
}

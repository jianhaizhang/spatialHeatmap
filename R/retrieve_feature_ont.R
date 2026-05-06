#' Retrieve spatial feature ontology for svgs
#'
#' @param term_id An ontology id such as 'UBERON_0000955'.

#' @return A list ontology information.

#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' Ooms J (2014). "The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects." _arXiv:1403.2805 [stat.CO]_. <https://arxiv.org/abs/1403.2805>.

#' @importFrom utils URLencode

ols_get_term <- function(term_id) {
  
  stopifnot(is.character(term_id), length(term_id) == 1)
  
  # Convert OBO ID to IRI if needed
  if (!grepl("^http", term_id)) {
    term_id <- paste0(
      "http://purl.obolibrary.org/obo/",
      gsub(":", "_", term_id)
    )
  }
  
  base_url <- "https://www.ebi.ac.uk/ols4/api/terms?iri="
  url <- paste0(base_url, URLencode(term_id, reserved = TRUE))
  
  txt <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) return(NULL)
  )
  
  if (is.null(txt) || length(txt) == 0) return(NULL)
  
  pkg <- check_pkg('jsonlite'); if (is(pkg, 'character')) stop(pkg)
  json <- tryCatch(
    jsonlite::fromJSON(paste(txt, collapse = ""), simplifyVector = FALSE),
    error = function(e) return(NULL)
  )
  
  if (is.null(json) || is.null(json$`_embedded`) || is.null(json$`_embedded`$terms)) {
    return(NULL)
  }
  
  terms <- json$`_embedded`$terms
  if (length(terms) == 0) return(NULL)
  
  t <- terms[[1]]
  
  list(
    id = if (!is.null(t$obo_id)) t$obo_id else NA_character_,
    iri = if (!is.null(t$iri)) t$iri else NA_character_,
    label = if (!is.null(t$label)) t$label else NA_character_,
    description = if (!is.null(t$description)) {
      paste(t$description, collapse = " ")
    } else NA_character_,
    ontology = if (!is.null(t$ontology_name)) t$ontology_name else NA_character_
  )
}

#' Formatting an ontology ID
#'
#' @param term_id An ontology id such as 'UBERON_0000955'.

#' @return An ontology id in the form a URL.

#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

return_ols <- function(term_id) {
  
  # If OBO ID like UBERON:0000955 → convert to IRI
  if (!grepl("^http", term_id)) {
    term_id <- paste0(
      "http://purl.obolibrary.org/obo/",
      gsub(":", "_", term_id)
    )
  }
  
  ols_get_term(term_id)
}

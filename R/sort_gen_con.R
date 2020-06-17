#' Sort "gene_condition" Names  
#'
#' @param ID.sel The vector of target genes.
#' @param na.all All "gene_condition" names.
#' @param con.all All condition names.
#' @param by "gene" or "con". Sort conditions for each gene or vice versa. 

#' @return A vector of sorted "gene_condition" names.
#' @keywords Internal

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

sort_gen_con <- function(ID.sel, na.all, con.all, by='gene') {

  # Sort vector of letter and number mixture.
  sort_mix <- function(vec) {

    w <- is.na(as.numeric(gsub('\\D', '', vec)))
    let <- vec[w]; num <- vec[!w]
    let.num <- as.numeric(gsub('\\D', '', vec))
    vec[!w] <- num[order(let.num[!w])]
    vec[w] <- let[order(let.num[w])]; return(vec)
        
  }  

  if (by=='gene') {

    # Sort conditions under each gene.
    con.pat1 <- paste0('.*_(', paste0(con.all, collapse='|'), ')$')
    na.sort <- NULL; for (i in sort_mix(ID.sel)) {
        
      na0 <- na.all[grepl(paste0('^', i, '_'), na.all)]
      con1 <- gsub(con.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(i, '_', sort_mix(con1)))

    }

  } else if (by=='con') {

    # Sort conditions and genes.
    gen.pat1 <- paste0('^(', paste0(ID.sel, collapse='|'), ')_.*')
    na.sort <- NULL; for (i in sort_mix(con.all)) {
      
      na0 <- na.all[grepl(paste0('_', i, '$'), na.all)]
      gen1 <- gsub(gen.pat1, '\\1', na0)
      na.sort <- c(na.sort, paste0(sort_mix(gen1), '_', i))

    }

  }; return(na.sort)

}

#' Extract up and down DEGs from any pair of comparisons in the provided data frame.
#'
#' @param sam.all All the samples/spatial features compared with each other.
#' @param df.all The data frame of of all pairwise comparisons from edgeR, limma, or DESeq2.
#' @param log.fc The cutoff of log fold change.
#' @param fdr The cutoff of the FDR.
#' @param log.na The strings in a column name indicating log fold change such as "log2FoldChange".
#' @param fdr.na The strings in a column name indicating FDR such as "padj".

#' @return A nested list of separated up and down genes.

#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

up_dn <- function(sam.all, df.all, log.fc, fdr, log.na, fdr.na, outliers=0, method=NULL) {
  # save(sam.all, df.all, log.fc, fdr, log.na, fdr.na, file='up.dn.arg')
  total <- NULL
  lis <- NULL; for (i in sam.all) {
    df.up <- df.up1 <- df.down <- df.down1 <- data.frame()
    w.fc <- grep(paste0('^', i, '_VS_.*_', log.na , '$'), colnames(df.all))
    w.fdr <- grep(paste0('^', i, '_VS_.*_', fdr.na, '$'), colnames(df.all))
    w.fc1 <- grep(paste0('_VS_', i, '_', log.na, '$'), colnames(df.all))
    w.fdr1 <- grep(paste0('_VS_', i, '_', fdr.na, '$'), colnames(df.all))
    if (length(w.fc)>0) {
      df.fc <- df.all[, w.fc, drop=FALSE]
      df.fdr <- df.all[, w.fdr, drop=FALSE]
      # Up DEGs.
      df.up.idx <- (df.fc >= abs(log.fc)) + (df.fdr <= fdr)
      # Down DEGs
      df.down.idx <- (df.fc <= -abs(log.fc)) + (df.fdr <= fdr)
    } # "else if" runs only if the preceeding "if" is not true (exclusive relationship), so should not be used here.

    if (length(w.fc1)>0) {
      df.fc <- df.all[, w.fc1, drop=FALSE]
      df.fdr <- df.all[, w.fdr1, drop=FALSE]
      # Up DEGs. 
      df.up.idx1 <- (df.fc <= -abs(log.fc)) + (df.fdr <= fdr) 
      # Down DEGs 
      df.down.idx1 <- (df.fc >= abs(log.fc)) + (df.fdr <= fdr)
    }

    if (length(w.fc)>0 & length(w.fc1)>0) {
      up.idx.all <- cbind(df.up.idx, df.up.idx1)
      dn.idx.all <- cbind(df.down.idx, df.down.idx1)
    } else if (length(w.fc)>0 & length(w.fc1)==0) {
      up.idx.all <- df.up.idx; dn.idx.all <- df.down.idx
    } else if (length(w.fc)==0 & length(w.fc1)>0) {
      up.idx.all <- df.up.idx1; dn.idx.all <- df.down.idx1
    }
    
    # Subset ups.
    rsum.up <- rowSums(up.idx.all==2)
    w.up <- which(rsum.up >= ncol(up.idx.all)-outliers)
    up <- df.all[w.up, c(w.fc, w.fc1, w.fdr, w.fdr1), drop=FALSE]
    # Include total references that ups are relative to.
    up <- cbind(total=rsum.up[w.up], up)

    # Subset dns.
    rsum.dn <- rowSums(dn.idx.all==2)       
    w.dn <- which(rsum.dn >= ncol(dn.idx.all)-outliers)
    dn <- df.all[w.dn, c(w.fc, w.fc1, w.fdr, w.fdr1), drop=FALSE]
    # Include total references that dns are relative to.
    dn <- cbind(total=rsum.dn[w.dn], dn)
    len <- length(sam.all)
    if (nrow(up)==0) up <- data.frame() else {
      up <- as.data.frame(up); cna.up <- colnames(up)
      fdr.up <- up[, grep(paste0('_', fdr.na, '$'), cna.up), drop=FALSE]
      up <- up[, order(cna.up)]
      up <- cbind(FDR_mean=10^rowMeans(log10(fdr.up)), type='up', method=method, up)
      cna.up <- colnames(up) # Necessary.
      cna.sel.up <- c('type', 'total', 'method', 'FDR_mean')
      up <- up[order(up[, 'FDR_mean']), c(cna.sel.up, setdiff(cna.up, cna.sel.up))]
      up <- rbind(subset(up, total==(len-1)), subset(up, total!=(len-1)))
    }

    if (nrow(dn)==0) dn <- data.frame() else {
      dn <- as.data.frame(dn); cna.dn <- colnames(dn)
      fdr.dn <- dn[, grep(paste0('_', fdr.na, '$'), cna.dn), drop=FALSE]
      dn <- dn[, order(colnames(dn))]
      dn <- cbind(FDR_mean=10^rowMeans(log10(fdr.dn)), type='down', method=method, dn)
      cna.dn <- colnames(dn) # Necessary.
      cna.sel.dn <- c('type', 'total', 'method', 'FDR_mean') 
      dn <- dn[order(dn[, 'FDR_mean']), c(cna.sel.dn, setdiff(cna.dn, cna.sel.dn))]
      dn <- rbind(subset(dn, total==(len-1)), subset(dn, total!=(len-1)))
    }; cat(i, 'up:', nrow(up), ';', 'down:', nrow(dn), '\n')

    lis0 <- list(up=up, down=dn)
    # names(lis0) <- paste0(i, c('_up', '_down'))
    lis1 <- list(lis0); names(lis1) <- i; lis <- c(lis, lis1)
   }; return(lis)
}

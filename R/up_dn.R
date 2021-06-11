#' Extract up and down DEGs from any pair of comparisons in the provided data frame.
#'
#' @param sam.all All the samples/spatial features compared with each other.
#' @param df.all The data frame of of all pairwise comparisons from edgeR, limma, DESeq2, or distinct.
#' @param log.fc The cutoff of log fold change.
#' @param fdr The cutoff of the FDR.
#' @param log.na The strings in a column name indicating log fold change such as "log2FoldChange".
#' @param fdr.na The strings in a column name indicating FDR such as "padj".

#' @return A nested list of separated up and down genes.

#' @keywords Internal
#' @noRd

#' @author Jianhai Zhang \email{jzhan067@@ucr.edu; zhang.jianhai@@hotmail.com} \cr Dr. Thomas Girke \email{thomas.girke@@ucr.edu}

#' @references 
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

up_dn <- function(sam.all, df.all, log.fc, fdr, log.na, fdr.na) {

  lis <- NULL; for (i in sam.all) {

    df.up <- df.up1 <- df.down <- df.down1 <- data.frame()
    w.fc <- grep(paste0('^', i, '_VS_.*_', log.na , '$'), colnames(df.all))
    w.fc1 <- grep(paste0('_VS_', i, '_', log.na, '$'), colnames(df.all))
    if (length(w.fc)>0) {

      df.fc <- df.all[, w.fc, drop=FALSE]
      w.fdr <- grep(paste0('^', i, '_VS_.*_', fdr.na, '$'), colnames(df.all))
      df.fdr <- df.all[, w.fdr, drop=FALSE]
      # Up DEGs.
      df.idx <- cbind(df.fc >= abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.up <- cbind(df.fc[w, , drop=FALSE], df.fdr[w, , drop=FALSE])
      }
      # Down DEGs
      df.idx <- cbind(df.fc <= -abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.down <- cbind(df.fc[w, , drop=FALSE], df.fdr[w, , drop=FALSE])
      }

    } # "else if" runs only if the preceeding "if" is not true (exclusive relationship), so should not be used here.

    if (length(w.fc1)>0) {

      df.fc <- df.all[, w.fc1, drop=FALSE]
      w.fdr1 <- grep(paste0('_VS_', i, '_', fdr.na, '$'), colnames(df.all))
      df.fdr <- df.all[, w.fdr1, drop=FALSE]
      # Up DEGs.
      df.idx <- cbind(df.fc <= -abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.up1 <- cbind(df.fc[w, , drop=FALSE], df.fdr[w, , drop=FALSE])
      }
      # Down DEGs
      df.idx <- cbind(df.fc >= abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.down1 <- cbind(df.fc[w, , drop=FALSE], df.fdr[w, , drop=FALSE])
      }

    }

    if (length(w.fc)>0 & length(w.fc1)>0) {

      up.int <- intersect(rownames(df.up), rownames(df.up1))
      down.int <- intersect(rownames(df.down), rownames(df.down1))
      up <- cbind(df.up[up.int, , drop=FALSE], df.up1[up.int, , drop=FALSE])
      down <- cbind(df.down[down.int, , drop=FALSE], df.down1[down.int, , drop=FALSE])
    } else { up <- rbind(df.up, df.up1); down <- rbind(df.down, df.down1) }

    if (nrow(up)==0) up <- data.frame() else {
      up <- as.data.frame(up)
      fdr.up <- up[, grep(paste0('_', fdr.na, '$'), colnames(up)), drop=FALSE]; up <- up[, order(colnames(up))]
      up <- cbind(FDR_mean=10^rowMeans(log10(fdr.up)), up); up <- up[order(up[, 1]), ]
    }

    if (nrow(down)==0) down <- data.frame() else {
      down <- as.data.frame(down)
      fdr.down <- down[, grep(paste0('_', fdr.na, '$'), colnames(down)), drop=FALSE]; down <- down[, order(colnames(down))]
      down <- cbind(FDR_mean=10^rowMeans(log10(fdr.down)), down); down <- down[order(down[, 1]), ]
    }; cat(i, 'up:', nrow(up), ';', 'down:', nrow(down), '\n')
    lis0 <- list(up=up, down=down); names(lis0) <- paste0(i, c('_up', '_down')); lis1 <- list(lis0); names(lis1) <- i
    lis <- c(lis, lis1)

   }; return(lis)

}

# Create FC matrix from real data. For a certain gene, the FC in one tissue >= de.min or <= 1/de.min is selected. If two FCs from two tissues fulfil reach criteria, the larger one between FC and 1/FC is selected.
fc_mat <- function(se, com.factor, de.min=2) {
  # All FCs.
  fc.all <- fc(se, com.factor=com.factor)
  # A vector all of FCs = 1.
  fc.sel <- rep(1, nrow(fc.all))
  # All max FCs and all min FCs.
  max.all <- vapply(seq_len(nrow(fc.all)), function(x) {max(fc.all[x, ])}, numeric(1))
  min.all <- vapply(seq_len(nrow(fc.all)), function(x) {min(fc.all[x, ])}, numeric(1))
  # Select FCs that is the only one >= de.min for a gene.
  fc.idx1 <- rowSums(fc.all >= de.min) == 1
  fc.sel[fc.idx1] <- max.all[fc.idx1]

  # Select FCs that is the only one <= 1/de.min for a gene.
  fc.idx2 <- rowSums(fc.all < 1/de.min) == 1
  fc.sel[fc.idx2] <- min.all[fc.idx2]
  # Index of genes whos FCs >= de.min and FCs <= 1/de.min.
  fc.idx3 <- fc.idx1 & fc.idx2
  # Compare FCs and 1/FCs and select the larger one.
  fc.max0 <- max.all[fc.idx3]; fc.min0 <- min.all[fc.idx3]
  fc.sel[fc.idx3] <- ifelse(fc.max0 > 1/fc.min0, fc.max0, fc.min0)
  # All FCs not selected are assigned 1.
  fc.all[(fc.all - fc.sel) != 0] <- 1.00
  fcts <- colData(se)[, com.factor]; fct <- fcts[!duplicated(fcts)]
  mat.fc <- fc.all[, rep(seq_len(ncol(fc.all)), table(fcts)[fct])]; return(mat.fc)
}

# Create FC matrix by randomly selecting FCs in each tissue, regardless if the selected FC is max or min in the corresponding gene.
fc_mat1 <- function(se, com.factor, de.p=0.2, fc.all, de.min=2, seed=1) {
  fcts <- colData(se)[, com.factor]; fct.uni <- fcts[!duplicated(fcts)]
  if (!all(sort(fct.uni) == sort(colnames(fc.all)))) stop('Make sure unique factors in "com.factor" are the same with column names in "fc.all"!')
  # Compute DEG ratios across tissues.
  set.seed(seed); mat.fc <- sam_fc(fc.all, fc.min=de.min, n=1)
  # Function to sample FCs for a single tissue.
  fc_mat0 <- function(fct) {
    # Number of DEGs.
    n <- round(nrow(se)*de.p*mat.fc['ratio', fct])
    # Exclude FCs equal to 0.
    fcs <- fc.all[, fct]; fcs <- fcs[fcs != 0]
    # Select FC >= provided DE.
    idx1 <- fcs >= 1; fcs1 <- fcs[idx1]
    idx2 <- fcs < 1; fcs2 <- fcs[idx2]
    fcs4 <- c(fcs1[fcs1 >= de.min], fcs2[fcs2 < 1/de.min])
    mat <- matrix(rep(1, n*ncol(fc.all)), nrow=n)
    # Sample FCs and create FC matrix for a single tissue.
    set.seed(seed); mat[, fct.uni == fct] <- sample(fcs4, n)
    return(mat)
  }
  # Combine all FC matrix.
  fc.mat <- do.call('rbind', lapply(fct.uni, fc_mat0))
  # The non-DEG matrix.
  mat.null <- matrix(rep(1, (nrow(se)-nrow(fc.mat))*3), ncol=3)
  mat.all <- rbind(fc.mat, mat.null)
  rownames(mat.all) <- rownames(se); colnames(mat.all) <- fct.uni
  # Complete FC matrix with replicates.
  return(mat.all[, rep(seq_len(ncol(mat.all)), table(fcts)[fct.uni])])
}


# Data matrix normalisation and plotting. 

library(edgeR); library(DESeq2)

library(reshape2); library(ggplot2)
plot_gen <- function(se, id, col.bg='lightblue', col.tar='purple', size=15, angle=45) { 
  
  mat <- assay(se)
  # convert to long format
  mat1 <- mat[c(setdiff(rownames(mat), id), id), ]; df.long <- reshape2::melt(mat1)
  cols <- rep(col.bg, nrow(mat1)); cols[rownames(mat1) %in% id] <- col.tar
  # The levels in melted data frame has the same order (left to right) with row order (top to bottom) in original data frame before melted.       
  # Colours map to variables in original data frame before melted.                                                                 
  # Possible: the colour order (left to right) matches with the row order (top to bottom) in original data frame before melted, but the coloured lined is plotted in the order of levels (left to right) in melted data frame.   
  g <- ggplot(data=df.long, aes(x=Var2, y=value, colour=Var1, group=Var1))+geom_line()+theme(legend.position="none")+              
scale_color_manual(values=cols)+labs(title="", x="Sample", y="Value")+theme(axis.text=element_text(size=size), axis.title=element_text(size=size, face="bold"), axis.text.x=element_text(angle=angle, hjust=1))
  if (any(id %in% rownames(mat1))) {
  
    df.long1 <- subset(df.long, Var1 %in% id) 
    g <- g+geom_point(data=df.long1,aes(x=Var2, y=value))+geom_text(data=df.long1, aes(x=Var2, y=value, label=round(value, 2), vjust=-0.5))+geom_text(data=subset(df.long1, Var2 %in% colnames(mat1)[ncol(mat1)]), aes(label=Var1, y=value, vjust=1, hjust=-0.2), angle=90)
  
  }; return(g)
 
}            

# Extract up and down DEGs from any pair of comparisons. If no argument of "df.all", function up_dn will look for 'df.all' from R console (global environment), not from within edgeR or deseq2. Therefore, it is necessary to set the "df.all" argument.
# Allow df.all to be a matrix or data frame.
up_dn <- function(sam.all, df.all=df.all, log.fc, fdr, log.na, fdr.na) {

  lis <- NULL; for (i in sam.all) {

    df.up <- df.up1 <- df.down <- df.down1 <- data.frame()
    w.fc <- grep(paste0('^', i, '_VS_.*_', log.na , '$'), colnames(df.all))
    w.fc1 <- grep(paste0('_VS_', i, '_', log.na, '$'), colnames(df.all))
    if (length(w.fc)>0) {

      df.fc <- df.all[, w.fc, drop=F]
      w.fdr <- grep(paste0('^', i, '_VS_.*_', fdr.na, '$'), colnames(df.all))
      df.fdr <- df.all[, w.fdr, drop=F]
      # Up DEGs.
      df.idx <- cbind(df.fc >= abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.up <- cbind(df.fc[w, , drop=F], df.fdr[w, , drop=F])
      }
      # Down DEGs
      df.idx <- cbind(df.fc <= -abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.down <- cbind(df.fc[w, , drop=F], df.fdr[w, , drop=F])
      }

    } # "else if" runs only if the preceeding "if" is not true (exclusive relationship), so should not be used here.

    if (length(w.fc1)>0) {

      df.fc <- df.all[, w.fc1, drop=F]
      w.fdr1 <- grep(paste0('_VS_', i, '_', fdr.na, '$'), colnames(df.all))
      df.fdr <- df.all[, w.fdr1, drop=F]
      # Up DEGs.
      df.idx <- cbind(df.fc <= -abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.up1 <- cbind(df.fc[w, , drop=F], df.fdr[w, , drop=F])
      }
      # Down DEGs
      df.idx <- cbind(df.fc >= abs(log.fc), df.fdr <= fdr)
      w <- which(rowSums(df.idx)==ncol(df.idx))
      if (length(w)>0) {
        df.down1 <- cbind(df.fc[w, , drop=F], df.fdr[w, , drop=F])
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
      fdr.up <- up[, grep(paste0('_', fdr.na, '$'), colnames(up)), drop=F]; up <- up[, order(colnames(up))]
      up <- cbind(FDR_mean=10^rowMeans(log10(fdr.up)), up); up <- up[order(up[, 1]), ]
    }

    if (nrow(down)==0) down <- data.frame() else {
      down <- as.data.frame(down)
      fdr.down <- down[, grep(paste0('_', fdr.na, '$'), colnames(down)), drop=F]; down <- down[, order(colnames(down))]
      down <- cbind(FDR_mean=10^rowMeans(log10(fdr.down)), down); down <- down[order(down[, 1]), ]
    }; cat(i, 'up:', nrow(up), ';', 'down:', nrow(down), '\n')
    lis0 <- list(up=up, down=down); names(lis0) <- paste0(i, c('_up', '_down')); lis1 <- list(lis0); names(lis1) <- i
    lis <- c(lis, lis1)

   }; return(lis)

}

# log2FC of con1 VS con2 for gene1 in the output approximates to log2(mean(con1)/mean(con2)) in the un-normalized count matrix for gene1.
library(edgeR)
edgeR <- function(se, method.norm='TMM', com.factor, method.adjust='BH', return.all=FALSE, log2.fc=1, fdr=0.05) {

  df.cnt <- assay(se); y <- DGEList(counts=df.cnt)
  cat('Normalizing ...', '\n') # norm.factors are used in glmFit.
  # To store normalized counts in 'se' and set method.norm='none' is not right, since the 'norm.factors' are essentially used but they are 1 if method.norm='none'.
  y <- calcNormFactors(y, method=method.norm)
  fct <- factor(colData(se)[, com.factor]); design <- model.matrix(~0+fct)
  colnames(design) <- levels(fct); rownames(design) <- colnames(df.cnt) # Duplicated row names in design are allowed.
  y <- estimateDisp(y, design); fit <- glmFit(y, design)

  com <- combn(x=colnames(design), m=2); con <- paste(com[1,], com[2,], sep="-")
  con.ma <- makeContrasts(contrasts=con, levels=design)
  cna.con <- colnames(con.ma); cna.con1 <- sub('-', '_VS_', cna.con)
  cat('Computing DEGs ...', '\n')
  df.all <- data.frame(rm=rep(NA, nrow(df.cnt))); for (i in seq_len(length(cna.con))) {

    lrt <- glmLRT(fit, contrast=con.ma[, i])
    tag <- as.data.frame(topTags(lrt, n=nrow(lrt), adjust.method=method.adjust, sort.by="none", p.value=1))
    cat(cna.con1[i], '\n')
    colnames(tag) <- paste0(cna.con1[i], '_', colnames(tag)); df.all <- cbind(df.all, tag)

  }; df.all <- df.all[, -1]
  if (return.all==TRUE) return(df.all)
  UD <- up_dn(sam.all=levels(fct), df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='logFC', fdr.na='FDR'); return(UD)

}

# log2FC of con1 VS con2 for gene1 in the output approximates to log2(mean(con1)/mean(con2)) in the un-normalized count matrix for gene1.
library(DESeq2)
deseq2 <- function(se, com.factor, method.adjust='BH', return.all=FALSE, log2.fc=1, fdr=0.05) {

  expr <- assay(se)
  col.data <- as.data.frame(colData(se)); fct <- col.data[, com.factor] <- as.factor(col.data[, com.factor])
  # Change the column name of target comparison elements to a uniform name, since user-provided names are different.   
  colnames(col.data)[colnames(col.data)==com.factor] <- 'Factor'

  # Alternatively, if use a design matrix for design, then use a numeric contrast vector to contrasts in "results", such as c(-1,  0, 0, 1). DESeq2 operates on raw counts.
  dds <- DESeqDataSetFromMatrix(countData=expr, colData=col.data, design=~Factor) # "design" does not impact "rlog" and "varianceStabilizingTransformation". Rownames of colData do not need to be identical with countData.
  dds <- DESeq(dds)

  com <- combn(x=levels(fct), m=2)
  df.all <- data.frame(rm=rep(NA, nrow(expr))); for (i in seq_len(ncol(com))) {

    contr <- c('Factor', com[2, i], com[1, i]); cat(contr[2], '-', contr[3], '\n')
    df0 <- as.data.frame(results(object=dds, contrast=contr, pAdjustMethod=method.adjust)[, c('log2FoldChange', 'padj')])
    colnames(df0) <- paste0(contr[2], '_VS_', contr[3], '_', colnames(df0)); df.all <- cbind(df.all, df0)

  }; df.all <- df.all[, -1]
  if (return.all==TRUE) return(df.all)
  UD <- up_dn(sam.all=levels(fct), df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='log2FoldChange', fdr.na='padj'); return(UD)

}

# log2FC of con1 VS con2 for gene1 in the output approximates to log2(mean(con1)/mean(con2)) in the un-normalized count matrix for gene1.
library(limma)
limma <- function(se, m.array, method.norm='TMM', com.factor, method.adjust='BH', return.all=FALSE, log2.fc=1, fdr=0.05) {

  # Design matrix.
  expr <- assay(se); fct <- factor(colData(se)[, com.factor]); design <- model.matrix(~0+fct)
  colnames(design) <- levels(fct); rownames(design) <- colnames(expr)

  # RNA-seq data.
  if (m.array==FALSE) {
    # The limma tutorial recommends TMM normalization.
    y <- DGEList(counts=assay(se)); cat('Normalising:', method.norm, '\n')
    # To store normalized counts in 'se' and set method.norm='none' is not right, since the 'norm.factors' are essentially used but they are 1 if method.norm='none'.
    y <- calcNormFactors(y, method=method.norm)
    v <- voom(y, design, plot=FALSE); fit <- lmFit(v, design)

  } else if (m.array==TRUE) fit <- lmFit(expr, design) # Microarray data.

  # Contrast matrix.
  com <- combn(x=colnames(design), m=2); con <- paste(com[1,], com[2,], sep="-")
  con.mat <- makeContrasts(contrasts=con, levels=design)
  cna.con <- colnames(con.mat); cna.con1 <- sub('-', '_VS_', cna.con)
  fit1 <- contrasts.fit(fit, con.mat); fit2  <- eBayes(fit1)

  df.all <- data.frame(rm=rep(NA, nrow(expr))); for (i in seq_along(cna.con)) {

    top <- topTable(fit2, coef=cna.con[i], number=Inf, p.value=1, adjust.method="BH", lfc=0, sort.by='none')
    colnames(top) <- paste0(cna.con1[i], '_', colnames(top)); df.all <- cbind(df.all, top)

  }; df.all <- df.all[, -1]
  if (return.all==TRUE) return(df.all)
  # Up and down DEGs.
  UD <- up_dn(sam.all=levels(fct), df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='logFC', fdr.na='adj.P.Val'); return(UD)

}


# At least three samples should be compared, otherwise errors arise. 
library(TCC); library(spatialHeatmap)
roku <- function(data, norm.fun='CNF', parameter.list=NULL, com.factor, aggr='mean', n=1, log2fc=1, return.all=FALSE) {
  if (!is(data, 'SummarizedExperiment')) stop('The "data" should be a "SummarizedExperiment" object!')
  fct.cna <- as.vector(colData(data)[, com.factor])
  df0 <- assay(data); cnt <- all(ceiling(df0)==df0) & all(df0>=0)
  if (cnt==TRUE) {
    cat('Count matrix is log2-transformed ... \n')
    if (!is.null(norm.fun)) data <- norm_data(data, norm.fun, parameter.list, log2.trans=TRUE)
    cat('Differential expression based on log2 FC in normalised count matrix: \n')
    df1 <- assay(data); colnames(df1) <- fct.cna
    lis.de <- de_fc(2^df1, log2FC=log2fc, aggr.rep=aggr)
  } else lis.de <- NULL
  df.exp <- assay(data); colnames(df.exp) <- fct.cna
  
  if (any(duplicated(fct.cna))) { cat('Aggregating replicates... \n') 
    # To keep colnames, "X" should be a character, not a factor.
    if (aggr=='mean') df.exp <- vapply(unique(fct.cna), function(x) rowMeans(df.exp[, fct.cna==x, drop=FALSE]), numeric(nrow(df.exp)))
    if (aggr=='median') { # The data to aggregate needs to be "matrix" when using "rowMedians".
      df.exp <- vapply(unique(fct.cna), function(x) Biobase::rowMedians(df.exp[, fct.cna==x, drop=FALSE]), numeric(nrow(df.exp)))
      rownames(df.exp) <- rownames(data)
    }
  }
  min.v <- min(df.exp); if (min.v<0) df.exp <- df.exp-min.v # Microarray data.
 
  # The data matrix should be non-negative values. 'n+0.5': make sure n is possible. E.g. 1/3=0.3333, if upper.limit is 0.3333, then max number of outlier is 3*0.3333 < 1.
  res <- ROKU(df.exp, upper.limit=(n+0.5)/ncol(df.exp)); res1 <- res[['outlier']]; res1 <- as.data.frame(cbind(res1, 
modH=res[['modH']]))
  if (return.all==TRUE) return(res1) else return(up_dn_roku(res1, lis.de))
}

# Separate up/down genes from ROKU results.
up_dn_roku <- function(df.all, lis.de=NULL) {                                                                                                   
  cna.all <- colnames(df.all); sam.all <- cna.all[-length(cna.all)]
  cat('Differential expression by ROKU: \n')                                                               
  lis <- NULL; for (i in sam.all) { 
    up0 <- subset(df.all, eval(parse(text=i))==1)[, c(i, 'modH')]                                                                
    if (nrow(up0)>0) {
      up0 <- up0[order(up0[, 'modH']), ]
      if (!is.null(lis.de)) up0 <- subset(up0, rownames(up0) %in% lis.de$up)
    }
    dn0 <- subset(df.all, eval(parse(text=i))==-1)[, c(i, 'modH')]                                                                 
    if (nrow(dn0)>0) {
      dn0 <- dn0[order(dn0[, 'modH']), ]
      if (!is.null(lis.de)) dn0 <- subset(dn0, rownames(dn0) %in% lis.de$down)
    }
    cat(i, 'up:', nrow(up0), 'down:', nrow(dn0), '\n') # No need to select up/down genes by the modH.
    lis0 <- list(up0, dn0); names(lis0) <- paste0(i, c('_up', '_down'))
    lis1 <- list(lis0); names(lis1) <- i; lis <- c(lis, lis1)
  }; return(lis)                                                                                                                   
}

library(distinct)
# Requires three columns in colData: sample, replicates, condition. All condition pairs are compared for each sample.
# Compare samples across all conditions: con.factor=NULL, which means actual samples are treated as conditions internally, and the sample column is assumed to include only one sample. Conversely, compare conditions across all samples: sam.factor=NULL, which means actual samples are treated as one sample internally.
# log2FC of con1 VS con2 for gene1 in the output approximates to log2(mean(con1)/mean(con2)) in the un-normalized count matrix for gene1.
distt <- function(se, norm.fun='CNF', parameter.list=NULL, log2.trans=TRUE, com.factor, return.all=FALSE, log2.fc=1, fdr=0.05) {

  if (!is.null(norm.fun)) { cat('Normalizing counts ... \n')
  se <- norm_data(data=se, norm.fun=norm.fun, parameter.list=parameter.list, log2.trans=log2.trans) }
   
  # Replicate names are the rownames of design matrix, and should not be duplicated, so avoid duplicated replicate names. 
  # cdat <- as.data.frame(colData(se)) # Duplicated row names are numbered.
  cdat <- colData(se); cdat$OneSample <- 'sample'; sam.factor <- 'OneSample'
  con.factor <- com.factor

  if (any(duplicated(rownames(cdat)))) rownames(cdat) <- paste0(rownames(cdat), '.', seq_along(rownames(cdat)))
  cdat$rep.distt <- rownames(cdat); colData(se) <- DataFrame(cdat)
  sam.all <- unique(cdat[, sam.factor])
  # Compare each pair of conditions for each sample.
  lis.all <- list(); for (i in sam.all) {
    idx <- cdat[, sam.factor] %in% i
    cdat0 <- cdat[idx, , drop=FALSE]
    con.all0 <- cdat0[, con.factor]; con.uni0 <- unique(con.all0)
    # Two or more conditions are required.
    if (length(unique(con.all0))==1) { cat(i, 'is skipped, since only one condition is detected! \n'); next }
    com <- combn(con.uni0, 2); rna <- rownames(se) 
    df.all1 <- data.frame(rm=rep(NA, length(rna)), row.names=rna)
    # Compare each condition pair of a sample.
    for (j in seq_len(ncol(com))) {
      con.all1 <- com[, j]
      idx1 <- (cdat[, con.factor] %in% con.all1) & idx
      se0 <- se[, idx1]
      fct <- factor(colData(se0)[, con.factor]); dsg <- model.matrix(~fct)
      vs <- paste0(levels(fct), collapse='_VS_')
      cat(i, ':', vs, '... \n')
      rownames(dsg) <- colData(se0)[, 'rep.distt']
      set.seed(61217)
      # Genes with values <=0 in less than min_non_zero_cells are disgarded.
      res0 <- distinct_test(x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_sample='rep.distt', design=dsg, column_to_test=2, min_non_zero_cells=ceiling((min(table(fct))+1)/2), n_cores=2)
      # Log2 FC requires count matrix.
      assay(se0) <- 2^assay(se0)
      res1 <- log2_FC(res=res0, x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_group=con.factor)
      cna.res <- colnames(res1)
      colnames(res1)[cna.res=='p_adj.loc'] <- paste0(vs, '_FDR') 
      colnames(res1)[grepl('^log2FC', cna.res)] <- paste0(vs, '_log2FC')
      rownames(res1) <- res1$gene; res1 <- subset(res1, filtered==FALSE)
      res1 <- res1[, grepl('_FDR$|_log2FC$|^filtered$', colnames(res1))]
      # Avoid NA in rownames.
      int0 <- intersect(rownames(df.all1), rownames(res1))
      df.all1 <- cbind(df.all1[int0, ], res1[int0, ])
    }; df.all1 <- df.all1[, -1]
    if (return.all==TRUE) return(df.all1)
    UD <- up_dn(sam.all=con.uni0, df.all=df.all1, log.fc=abs(log2.fc), fdr=fdr, log.na='log2FC', fdr.na='FDR'); cat('\n')
    lis0 <- list(UD); names(lis0) <- i; lis.all <- c(lis.all, lis0)
  }; return(lis.all)
}


distt1 <- function(se, sam.factor, con.factor, rep, log2.fc=1, fdr=0.05) {

  # Replicate names are the rownames of design matrix, and should not be duplicated, so avoid duplicated replicate names. 
# cdat <- as.data.frame(colData(se)) # Duplicated row names are numbered.
  cdat <- colData(se); cdat$OneSample <- 'sample'
  if (is.null(con.factor)&!is.null(sam.factor)) { # Compare samples across different coditions.
    con.factor <- sam.factor; sam.factor <- 'OneSample' 
  } else if (!is.null(con.factor)&is.null(sam.factor)) {
    sam.factor <- 'OneSample'
  }
  if (any(duplicated(rownames(cdat)))) rownames(cdat) <- paste0(rownames(cdat), '.', seq_along(rownames(cdat)))
  cdat[, rep] <- rownames(cdat); colData(se) <- DataFrame(cdat)
  sam.all <- unique(cdat[, sam.factor])
  # Compare each pair of conditions for each sample.
  lis.all <- list(); for (i in sam.all) {
    idx <- cdat[, sam.factor] %in% i
    cdat0 <- cdat[idx, , drop=FALSE]
    con.all0 <- cdat0[, con.factor]; con.uni0 <- unique(con.all0)
    # Two or more conditions are required.
    if (length(unique(con.all0))==1) { cat(i, 'is skipped, since only one condition is detected! \n'); next }
    com <- combn(con.uni0, 2); rna <- rownames(se) 
    df.all1 <- data.frame(rm=rep(NA, length(rna)), row.names=rna)
    # Compare each condition pair of a sample.
    for (j in seq_len(ncol(com))) {
      con.all1 <- com[, j]
      idx1 <- (cdat[, con.factor] %in% con.all1) & idx
      se0 <- se[, idx1]
      fct <- factor(colData(se0)[, con.factor]); dsg <- model.matrix(~fct)
      vs <- paste0(levels(fct), collapse='_VS_')
      cat(i, ':', vs, '... \n')
      rownames(dsg) <- colData(se0)[, rep]
      set.seed(61217)
      # Genes with values <=0 in less than min_non_zero_cells are disgarded.
      res0 <- distinct_test(x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_sample=rep, design=dsg, column_to_test=2, min_non_zero_cells=ceiling((min(table(fct))+1)/2), n_cores=2)
      # Log2 FC requires count matrix.
      assay(se0) <- 2^assay(se0)
      res1 <- log2_FC(res=res0, x=se0, name_assays_expression=assayNames(se0)[1], name_cluster=sam.factor, name_group=con.factor)
      cna.res <- colnames(res1)
      colnames(res1)[cna.res=='p_adj.loc'] <- paste0(vs, '_FDR') 
      colnames(res1)[grepl('^log2FC', cna.res)] <- paste0(vs, '_log2FC')
      rownames(res1) <- res1$gene; res1 <- subset(res1, filtered==FALSE)
      res1 <- res1[, grepl('_FDR$|_log2FC$|^filtered$', colnames(res1))]
      # Avoid NA in rownames.
      int0 <- intersect(rownames(df.all1), rownames(res1))
      df.all1 <- cbind(df.all1[int0, ], res1[int0, ])
    }; df.all1 <- df.all1[, -1]
    UD <- up_dn(sam.all=con.uni0, df.all=df.all1, log.fc=abs(log2.fc), fdr=fdr, log.na='log2FC', fdr.na='FDR'); cat('\n')
    lis0 <- list(UD); names(lis0) <- i; lis.all <- c(lis.all, lis0)
  }; return(lis.all)
}


# Select DE genes according to log2FC on counts.
de_fc <- function(mat, log2FC=1, aggr.rep='mean') {
  fct.cna <- colnames(mat) 
  # To keep colnames, "X" should be a character, not a factor.
  if (aggr.rep=='mean') mat.aggr <- vapply(unique(fct.cna), function(x) rowMeans(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
  if (aggr.rep=='median') {
    mat.aggr <- vapply(unique(fct.cna), function(x) Biobase::rowMedians(mat[, fct.cna==x, drop=FALSE]), numeric(nrow(mat)))
    rownames(mat.aggr) <- rownames(mat)
  }
  # Recombination of column factors.
  com <- combn(unique(fct.cna), 2)
  df.fc <- data.frame(rm=rep(NA, nrow(mat.aggr)))
  cna.fc <- cna.p <- NULL
  # All pairs of column factors.
  for (i in seq_len(ncol(com))) {
    s1 <- com[, i][1]; s2 <- com[, i][2]
    fc0 <- log2(mat.aggr[, s1]/mat.aggr[, s2])
    cna.fc <- c(cna.fc, paste0(s1, '_VS_', s2, '_log2FC'))
    cna.p <- c(cna.p, paste0(s1, '_VS_', s2, '_Pvalue'))
    df.fc <- cbind(df.fc, fc0)
  }
  df.fc <- df.fc[, -1, drop=FALSE]; colnames(df.fc) <- cna.fc
  df.p <- matrix(rep(1, nrow(df.fc)*ncol(df.fc)), nrow=nrow(df.fc))
  colnames(df.p) <- cna.p; df.vs <- cbind(df.fc, df.p)
  # Select genes with aggregated values above 'log2FC'.
  UD <- up_dn(sam.all=unique(fct.cna), df.all=df.vs, log.fc=abs(log2FC), fdr=1, log.na='log2FC', fdr.na='Pvalue')
  ups <- dns <- NULL; for (i in seq_along(UD)) {
    lis0 <- UD[[i]]; ups <- c(ups, rownames(lis0[[1]]))
    dns <- c(dns, rownames(lis0[[2]]))
  }; ups <- unique(ups); dns <- unique(dns); uni <- unique(c(ups, dns))
  # One gene can be up in sampleA while down in sampleB, so ups and dns can overlap.
  cat('Up:', length(ups), 'Down:', length(dns), 'Unique:', length(uni), '\n')
  return(list(up=ups, down=dns, unique=uni))
}

# FC of one condition or sample VS all other conditions or samples respectively. Inf/-Inf is avoided by adding 1 to all counts if the min count is 0 in the count matrix.
fc <- function(data, norm.fun='none', parameter.list=NULL, com.factor, aggr='mean') {
  if (!is(data, 'SummarizedExperiment')) stop('The "data" should be a "SummarizedExperiment" object!')
  fct.cna <- as.vector(colData(data)[, com.factor])
  df0 <- assay(data)
  if (!all(ceiling(df0)==df0) & all(df0>=0)) stop('The count matrix should only include non-negative values!')
  if (norm.fun!='none') { cat('Normalizing ... \n')
    data <- norm_data(data, norm.fun, parameter.list, log2.trans=FALSE)
  }
    mat <- assay(data); colnames(mat) <- fct.cna
    if (min(mat) == 0) mat <- mat + 1 # Avoid Inf values.
    if (aggr=='mean') mat.aggr <- vapply(unique(fct.cna), function(x) { 
      mn1 <- rowMeans(mat[, fct.cna==x, drop=FALSE])
      mn2 <- rowMeans(mat[, fct.cna!=x, drop=FALSE])
      return(mn1/mn2)
    }, numeric(nrow(mat)))

    if (aggr=='median') mat.aggr <- vapply(unique(fct.cna), function(x) { 
      md1 <- Biobase::rowMedians(mat[, fct.cna==x, drop=FALSE])
      md2 <- Biobase::rowMedians(mat[, fct.cna!=x, drop=FALSE])
      return(md1/md2)
    }, numeric(nrow(mat)))
    rownames(mat.aggr) <- rownames(data); return(mat.aggr)
}

# Randomly select FCs for each gene. 1/fc.min is also internally considered.
sam_fc <- function(fc.all, fc.min=2, n) {
  mat <- matrix(rep(0, ncol(fc.all)*n), nrow=n)
  rownames(mat) <- paste0('FC', seq_len(n))
  len <- NULL
  # All FCs in each sample/condition.
  for (i in seq_len(ncol(fc.all))) {
    # All FCs larger than fc.min and not Inf in each sample/condition.
    fc0 <- fc.all[, i]; fc0 <- fc0[c(abs(fc0) >= fc.min | abs(fc0) <= 1/fc.min) & abs(fc0)!=Inf]
    len <- c(len, length(fc0))
    # Randomly select n FCs for each sample/condition.
    mat[, i] <- sample(fc0, n, replace=FALSE)
  }; colnames(mat) <- colnames(fc.all)
  mat <- rbind(ratio=len/sum(len), mat)
  return(mat)
}

# Aggregate P values in the data frame of all pairwise comparison result.
p_aggr <- function(df.all, pat, aggr='median') {
  df.p <- df.all[, grepl(pat, colnames(df.all))]
  if (aggr=='median') p.aggr <- Biobase::rowMedians(as.matrix(df.p))
  if (aggr=='mean') p.aggr <- rowMeans(as.matrix(df.p))
  return(p.aggr)
}


# Given data and DEG data frame, return response-prediction pairs.

res_pre <- function(se, df.all, pat, aggr='median', response) {
  # In distinct, some genes are lost after DEG identification.
  left <- rownames(se) %in% rownames(df.all)
  se <- se[left, ]; response <- response[left]
  # Isolate and aggregate p values.
  df.p <- df.all[, grepl(pat, colnames(df.all))]
  if (aggr=='median') p.aggr <- Biobase::rowMedians(as.matrix(df.p))
  if (aggr=='mean') p.aggr <- rowMeans(as.matrix(df.p))
  # Valid genes.
  min.p <- min(p.aggr); p.aggr[p.aggr==0] <- min.p; idx <- !is.na(p.aggr)
  return(roc(response[idx], 1/p.aggr[idx]))
} 


# ROC of provided function on provided data.

roc_fun <- function(se.sim, lis.fun, label) {
  fun.na <- lis.fun[['fun']] # Function name.
  # Equivalent arguments of se, com.factor in the provided function. 
  lis0 <- list(se.sim, 'com.factor')
  names(lis0) <- c(lis.fun$se, lis.fun$com.factor)
  # Local arguments of the provided function.
  lis.fun1 <- lis.fun[!names(lis.fun) %in% c('fun', 'se', 'com.factor', 'pat')]
  # ROC
  res <- do.call(fun.na, c(lis0, lis.fun1))
  roc0 <- res_pre(se.sim, res, lis.fun[['pat']], 'median', label)
  # Coordinates of the ROCs, which will be used in ggplot2.
  df0 <- data.frame(threshold=roc0$threshold, specificity=roc0$specificities, sensitivity=roc0$sensitivities, method=fun.na)
  lis.f <- list(roc0, df0); names(lis.f) <- c(paste0('roc.', fun.na), 'coordinate')  
  return(lis.f)
}
   

# Given real data and a list of functions (including arguments), run simulation and performance evaluation iteratively.

com_perfm <- function(se, com.factor, n.round=10, n.gene, p.DEG=0.1, fc.min=2, ratio.fc, fun.lis, color, line.size=0.5) {
  fc.all <- round(fc(se, com.factor=com.factor), 1)
  if (missing(ratio.fc)) { mat.fc <- sam_fc(fc.all, fc.min=fc.min, n=n.round) } else { mat.fc <- ratio.fc; n.round <- nrow(mat.fc)-1 }
  g <- ggplot(mapping=aes(x=specificity, y=sensitivity, group=method))+scale_x_reverse()
  fun.na <- vapply(lis.funs, function(x) return(x$fun), character(1))
  g.lis <- rep(list(NA), length(fun.lis))
  for (x in seq_along(g.lis)) g.lis[[x]] <- g
  auc.lis <- rep(list(NULL), length(fun.lis))
  names(g.lis) <- names(auc.lis) <- fun.na

  if (missing(color)) {
    cols <- grDevices::colors()[grep('honeydew|aliceblue|white|gr(a|e)y', grDevices::colors(), invert=TRUE)]
    n <- length(fun.lis) 
    color <- cols[seq(from=1, to=length(cols), by=floor(length(cols)/n))][seq_along(fun.lis)]
  }
  for (i in seq_len(n.round)) {
    cat('Round', i, '... \n'); fcts <- colData(se)[, com.factor]
    fct <- fcts[!duplicated(fcts)]
    tcc <- simulateReadCounts(mat=assay(se), Ngene=ifelse(missing(n.gene), nrow(se), n.gene), PDEG=p.DEG, DEG.assign=mat.fc['ratio', ], DEG.foldchange=mat.fc[i+1, ], replicates=table(fcts)[fct])

    mat <- tcc$count
    cdat <- data.frame(com.factor=rep(paste0('group', seq_along(fct)), each=3))
    se.sim <- SummarizedExperiment(assays=list(count=as.matrix(mat)), colData=cdat)
    des <- tcc$simulation$trueDEG; des <- names(des)[des==1]

    rna <- rownames(se.sim); label <- rep(0, length(rna))
    label[rna %in% des] <- 1

    df.cor <- data.frame(); for (j in seq_along(fun.lis)) {
      roc.cor <- roc_fun(se.sim, fun.lis[[j]], label) 
      df0 <- roc.cor$coordinate
      df0$method <- paste0(df0$method, '.', i)
      g.lis[[j]] <- g.lis[[j]]+geom_line(data=df0, linetype="solid", color=color[j], size=line.size)+theme(legend.position=c(0.85, 0.1), plot.title=element_text(hjust=0.5))+ggtitle(paste0('ROCs of ', sub('\\.\\d+$', '', df0$method[1]))) 
      if (n==1) g.lis[[j]] <- g.lis[[j]]+theme(plot.title=element_text(hjust=0.5))+ggtitle(paste0('ROCs of ', sub('\\.\\d+$', '', df0$method[1]))) 
      auc.lis[[j]] <- c(auc.lis[[j]], as.numeric(auc(roc.cor[[1]])))
      df.cor <- rbind(df.cor, df0)
    }
    funs <- unique(df.cor$method); col.all <- rep(NA, nrow(df.cor))
    for (k in seq_along(funs)) col.all[df.cor$method %in% funs[k]] <- color[k]
    df.cor$method <- factor(df.cor$method, levels=unique(df.cor$method))
    # aes(color=var) is needed to make a legend, which means color variables according to elements in var, then use scale_color_manual(values=) to change colors.
    if (i==1) {
      names(col.all) <- df.cor$method
      g <- g+geom_line(data=df.cor, aes(color=method), linetype="solid", size=line.size)+scale_color_manual(values=col.all, labels=sub('\\.\\d+$', '', levels(df.cor$method)))+theme(legend.position=c(0.85, 0.1), plot.title=element_text(hjust=0.5))+ggtitle('ROCs of All Methods') } else {
      g <- g+geom_line(data=df.cor, linetype="solid", color=col.all, size=line.size)
     }

    df.auc <- data.frame(); for (i in seq_along(auc.lis)) {
      df0 <- cbind(method=names(auc.lis)[i], auc=round(auc.lis[[i]], 2))
      df.auc <- rbind(df.auc, df0)
    }
  } 
  df.auc$method <- factor(df.auc$method, levels=unique(df.auc$method))
  df.auc$auc <- as.numeric(df.auc$auc)
  return(list(roc.all=g, roc.each=g.lis, auc.all=df.auc))
}



# lis <- list(edg=deg.edg, dsq=deg.dsq, lim=deg.lim.r, tcc=deg.tcc.r)
# sam <- 'merzo.sb.os'
# Extract up or down DEGs for provided sample, which are fed to UpsetR.
deg_lis <- function(lis, sam, deg='up') {

  sam.all <- NULL; for (i in lis) { sam.all <- c(sam.all, names(i)) }
  sam.all <- unique(sam.all)
  meths <- names(lis); if (length(lis)!=length(lis)) stop('Each element in "lis" should have a name!')
  lis.all <- NULL; for (i in meths) {
    lis0 <- lis[[i]]
    rna <- rownames(lis0[[sam]][[paste0(sam, '_', deg)]])
    if (length(rna)>0) { lis1 <- list(rna); names(lis1) <- paste0(i, '.', deg); lis.all <- c(lis.all, lis1) } else cat('No genes detected in', i, '! \n')
  }; return(lis.all)

}


# library(RankProd)
rank_prod <- function(se, com.factor, logged, base, pfp=0.05, log2.fc=1) {

  # Use original expression values.
  if (logged==TRUE) df.exp <- base^assay(se) else df.exp <- assay(se)
  rna <- rownames(df.exp)
  fct <- as.character(colData(se)[, com.factor])
  sam.all <- unique(fct); com <- combn(sam.all, 2)
  df.all <- matrix(data=rep(NA, 2*nrow(df.exp)), ncol=1)
  vs1 <- vs2 <- NULL; for (i in seq_len(ncol(com))) {

    cl1 <- com[1, i]; cl2 <- com[2, i]; cat(cl1, '-', cl2, '\n')
    df1 <- df.exp[, fct==cl1]; df2 <- df.exp[, fct==cl2]
    df3 <- cbind(df1, df2)
    cl <- c(rep(0, ncol(df1)), rep(1, ncol(df2)))
    rp <- RankProducts(data=df3, cl=cl, logged=FALSE, na.rm=FALSE, gene.names=rna, rand=123, calculateProduct=TRUE)
    top <- topGene(rp, cutoff=NULL, method=NULL, logged=FALSE, logbase=NULL, gene.names=rna, num.gene=length(rna))

    # All input genes cl1 < cl2.
    top[[1]][, 'FC:(class1/class2)'] <- log2(top[[1]][, 'FC:(class1/class2)'])
    vs1 <- top[[1]][rna, ]; colnames(vs1) <- paste0(cl1, '_VS_', cl2, '_', make.names(colnames(vs1)))
    # All input genes cl1 > cl2.
    top[[2]][, 'FC:(class1/class2)'] <- log2(top[[2]][, 'FC:(class1/class2)'])
    vs2 <- top[[2]][rna, ]; colnames(vs2) <- paste0(cl1, '_VS_', cl2, '_', make.names(colnames(vs2)))

    # Make df.all as a matrix, since each gene is present twice in the row names. Row names must be same when cbind data frame or  matrix.
    df0 <- rbind(vs1, vs2); df.all <- cbind(df.all, df0)

  }; df.all <- df.all[, -1]

  UD <- up_dn(sam.all=sam.all, df.all=df.all, log.fc=abs(log2.fc), fdr=pfp, log.na='FC..class1.class2.', fdr.na='pfp'); return(UD)

}

# library(baySeq)
bay_seq <- function(se, com.factor, lib.scale='quantile', quant=0.75, sample.size=1e5, rep.aggr='mean', log2.fc=1, fdr=0.05, clust=NULL, ...) {
  
  df.exp <- assay(se) 
  fct <- as.character(colData(se)[, com.factor]) 
  sam.all <- unique(fct); com <- combn(sam.all, 2)
  df.all <- matrix(data=rep(NA, nrow(df.exp)), ncol=1) 
  for (i in seq_len(ncol(com))) { 
       
    cl1 <- com[1, i]; cl2 <- com[2, i]; cat(cl1, '-', cl2, '\n') 
    df1 <- df.exp[, fct==cl1]; df2 <- df.exp[, fct==cl2]
    df3 <- cbind(df1, df2); rna <- rownames(df3) 
    cD <- new("countData", data=df3, annotation=data.frame(gene=rownames(df3), stringsAsFactors=FALSE))                            
    # Additional parameters to be passed to the 'edgeR' calcNormFactors function: ...                                              
    libsizes(cD) <- getLibsizes(cD, estimationType=lib.scale, quantile=quant, ...)                                                 
    replicates(cD) <- reps <- as.factor(c(rep(cl1, ncol(df1)), rep(cl2, ncol(df2))))                                               
    NDE <- factor(rep(1, ncol(df3))); com.pair <- as.character(reps)                                                               
    groups(cD) <- list(NDE=NDE, com.pair=com.pair) 
    if(!is.null(clust)) { if(clust>=1) { require("parallel"); cl <- makeCluster(clust) } else cl <- NULL } else cl <- NULL
    cD <- getPriors.NB(cD, samplesize=sample.size, cl=cl) 
    cD <- getLikelihoods(cD, nullData=TRUE, cl=cl) 
    if(!is.null(cl)) stopCluster(cl) 
    top <- topCounts(cD, group='com.pair', decreasing=FALSE, normaliseData=TRUE, number=nrow(df3), FDR=1)  
   # Calculate logFC based on mean counts. 
    df4 <- top[, colnames(df3)] 
    df.avg <- aggregate(t(df4), by=list(as.character(reps)), FUN=rep.aggr); cna <- df.avg[, 1]; df.avg <- t(df.avg[, -1])          
    colnames(df.avg) <- cna; rownames(df.avg) <- rownames(top) <- top$gene; top <- cbind(logFC=log2(df.avg[, cl1]/df.avg[, cl2]),  
top, stringsAsFactors=FALSE) 
    top <- top[, !(colnames(top) %in% c('gene', colnames(df3)))]            
    colnames(top) <- paste0(cl1, "_VS_", cl2, '_', colnames(top))
    df.all <- cbind(df.all, top[rna, ], stringsAsFactors=FALSE)

 }; df.all <- df.all[, -1] 
 UD <- up_dn(sam.all=sam.all, df.all=df.all, log.fc=abs(log2.fc), fdr=fdr, log.na='logFC', fdr.na='FDR.com.pair'); return(UD)
                       
}     

# Summarise SSGs by within-method significance and cross-method frequency.
sig_frq <- function(lis.all, per=0.5, sam=NULL) {

 lis.sig <- lis.frq <- NULL; if (is.null(sam)) sam <- names(lis.all[[1]]) 
  # Summarise SSGs by within-method significance.
  for (j in sam) {

    up.g <- dn.g <- NULL; for (i in seq_along(lis.all)) {

      df.up <- lis.all[[i]][[j]][[1]]; rn <- round(nrow(df.up)*per)
      if (!is.null(df.up) & nrow(df.up) <= rn) rn <- nrow(df.up)
      if (nrow(df.up)>0) up.g <- c(up.g, rownames(df.up[1:rn, ]))
      df.dn <- lis.all[[i]][[j]][[2]]; rn <- round(nrow(df.dn)*per)
      if (!is.null(df.dn) & nrow(df.dn) <= rn) rn <- nrow(df.dn)
      if (nrow(df.dn)>0) dn.g <- c(dn.g, rownames(df.dn[1:rn, ]))

    }; up.g <- sort(table(up.g), decreasing=TRUE); dn.g <- sort(table(dn.g), decreasing=TRUE)

    lis <- list(up.g, dn.g); names(lis) <- paste0(j, c('_up_sig', '_down_sig'))
    lis1 <- list(lis); names(lis1) <- j; lis.sig <- c(lis.sig, lis1)

  }

  # Summarise SSGs by across-method frequency.
  for (j in sam) {

    up.g <- dn.g <- NULL; for (i in seq_along(lis.all)) {
 
      up.g <- c(up.g, rownames(lis.all[[i]][[j]][[1]]))
      dn.g <- c(dn.g, rownames(lis.all[[i]][[j]][[2]]))

    }; up.g <- sort(table(up.g), decreasing=TRUE); dn.g <- sort(table(dn.g), decreasing=TRUE)

    lis <- list(up.g, dn.g); names(lis) <- paste0(j, c('_up_frq', '_down_frq'))
    lis1 <- list(lis); names(lis1) <- j; lis.frq <- c(lis.frq, lis1)

  }; return(list(lis.sig=lis.sig, lis.frq=lis.frq))


}

# Plot SSGs separately. 
library(ggplot2)
ssg_sep <- function(lis.all, main=NULL, sam=NULL, by.sample=TRUE, angle=45, width=0.8, size=15, cols=NULL) { 

  if (is.null(sam)) sam.all <- names(lis.all[[1]]) else sam.all <- sam
  df.up.dn <- data.frame(); for (i in names(lis.all)) {

    for (j in sam.all) {

      up0 <- nrow(lis.all[[i]][[j]][[1]]); dn0 <- nrow(lis.all[[i]][[j]][[2]])
      df.up.dn <- rbind(df.up.dn, data.frame(sample=paste0(c(j, j), '_', c(i, i)), method=c('up', 'down'), count=c(up0, dn0)))
 
    }

  }; colnames(df.up.dn) <- c('sample', 'method', 'count')

  n <- length(unique(df.up.dn$method))
  if (is.null(cols)) {
 
    cols <- colors()[grep('honeydew|aliceblue|white|gr(a|e)y', colors(), invert=TRUE)]
    cols <- cols[seq(from=1, to=length(cols), by=floor(length(cols)/n))]
    
  }

  if (by.sample==TRUE) df.up.dn <- df.up.dn[order(df.up.dn$sample), ] 
  g <- ggplot(df.up.dn, aes(x=sample, y=count, fill=method))+geom_bar(stat="identity", position='stack', width=width)+scale_fill_manual(values=cols)+geom_text(position=position_stack(vjust=0.1), aes(y=count+mean(count)*0.05, label=count, hjust=0), angle=0)+theme(legend.position="right", axis.text=element_text(size=size), axis.title=element_text(size=size, face="bold"), axis.text.x=element_text(angle=angle, hjust=1), plot.title=element_text(hjust=0.5))+scale_x_discrete(limits=unique(df.up.dn$sample))+labs(title=main, x="sample", y="Number of genes"); return(g)

}

# Plot SSGs by within-method significance and cross-method frequency.
library(gridExtra)
ssg_sum <- function(lis.aggr, main.sig='Summarise SSGs by \n within-method significance', main.frq='Summarise SSGs by \n cross-method frequency', cols=NULL, width=0.85, size=15, angle=45) {

  df.aggr <- function (lis.aggr) {
    
    sam.all <- names(lis.aggr)
    df.frq <- data.frame(); for (i in sam.all) {

      df1 <- as.data.frame(table(lis.aggr[[i]][[1]]))
      df2 <- as.data.frame(table(lis.aggr[[i]][[2]]))
      if (nrow(df1)>0) {

        df1$Var1 <- paste0('up.', df1$Var1); df1 <- cbind(Sample=i, df1); df.frq <- rbind(df.frq, df1)

      }

      if (nrow(df2)>0) {

        df2$Var1 <- paste0('down.', df2$Var1); df2 <- cbind(Sample=i, df2); df.frq <- rbind(df.frq, df2)

      }

    }; return(df.frq)

  }

  plot_sum <- function(df.frq, size, angle, main, cols) {

    n <- length(unique(df.frq$Var1))
    # library(RColorBrewer)
    # qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    # col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))); cols <- col_vector[1:n]
    if (is.null(cols)) {
 
      cols <- colors()[grep('honeydew|aliceblue|white|gr(a|e)y', colors(), invert=TRUE)]
      cols <- cols[seq(from=1, to=length(cols), by=floor(length(cols)/n))]
    
    }

    g <- ggplot(df.frq, aes(x=Sample, y=Freq, fill=Var1))+geom_bar(stat="identity", position=position_dodge(width=width), width=width)+scale_fill_manual(values=cols)+geom_text(position=position_dodge(width=width), aes(y=Freq+mean(Freq)*0.05, label=Var1, hjust=0), angle=90)+theme(legend.position="none", axis.text=element_text(size=size), axis.title=element_text(size=size, face="bold"), axis.text.x=element_text(angle=angle, hjust=1), plot.title=element_text(hjust=0.5))+scale_x_discrete(limits=unique(df.frq$Sample))+labs(title=main, x="Sample", y="Number of genes")+coord_cartesian(ylim=c(0, max(df.frq[, 'Freq'])*1.1)); return(g)
    
    }
    df.frq1 <- df.aggr(lis.aggr=lis.aggr[['lis.sig']])
    g1 <- plot_sum(df.frq=df.frq1, size=size, angle=angle, main=main.sig, cols=cols)
    df.frq2 <- df.aggr(lis.aggr=lis.aggr[['lis.frq']])
    g2 <- plot_sum(df.frq=df.frq2, size=size, angle=angle, main=main.frq, cols=cols)    
    grid.arrange(g1, g2, nrow=1)

}

# TissueEnrich; HPA (Uhlén et al. 2015)
HPA <- function(se, fold=5, min.expr=1, rep.aggr=NULL, com.factor=NULL) {

  expData <- assay(se)
  if (!is.null(rep.aggr) & !is.null(com.factor)) {

    cat('Aggregating replicates', '\n')
    mat <- aggregate(x=t(expData), by=list(sam.var=colData(se)[, com.factor]), FUN=rep.aggr)
    sam.var <- mat[, 1];  mat <- t(mat[, -1]); colnames(mat) <- sam.var; expData <- mat

  }
  # Up.
  df <- c(); for (i in rownames(expData)) {

    tpm <- expData[i, ]; tpm <- sort(tpm, decreasing=TRUE); high <- tpm[1]
    if (high >= min.expr) {

      fc <- high/tpm[2]; if (is.na(fc)) next
      if (fc >= fold) { df <- rbind(df, c(i, names(tpm)[1], fc)) }

    }

  }; df <- as.data.frame(df); if (nrow(df)>0) { colnames(df) <- c('gene', 'tissue', 'fold.min'); df[, 3] <- as.numeric(df[, 3]) }
  # Down.
  df1 <- c(); for (i in rownames(expData)) {

    tpm <- expData[i, ]; tpm <- sort(tpm, decreasing=FALSE); low <- tpm[1]
    if (low >= min.expr) {

      fc <- tpm[2]/low; if (is.na(fc)) next
      if (fc >= fold) { df <- rbind(df1, c(i, names(tpm)[1], 1/fc)) }

    }

  }; df1 <- as.data.frame(df1); if (nrow(df1)>0) { colnames(df1) <- c('gene', 'tissue', 'fold.max'); df1[, 3] <- as.numeric(df1[, 3]) }

  tis <- unique(colnames(expData)); up0 <- down0 <- data.frame()
  lis <- NULL; for (i in tis) {

   if (nrow(df)>0) { 

     up0 <- subset(df, tissue==i)[, -2]; rownames(up0) <- up0[, 1]; up0 <- up0[, c(2, 1)]
     up0 <- up0[order(up0[, 1]), ]

   }

   if (nrow(df1)>0) {

     down0 <- subset(df1, tissue==i)[, -2]; rownames(down0) <- down0[, 1]; down0 <- down0[, c(2, 1)]
     down0 <- down0[order(down0[, 1]), ]

   }; cat(i, ' up: ', nrow(up0), '; ', 'down: ', nrow(down0), '\n')
   lis0 <- list(up=up0, down=down0); names(lis0) <- paste0(i, c('_up', '_down')); lis1 <- list(lis0); names(lis1) <- i
   lis <- c(lis, lis1)

  }; return(lis)

}










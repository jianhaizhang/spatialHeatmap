# ucr_efp_expr_ann_row_gen.txt

# Source: http://efp.ucr.edu/data/Mean_Intensities_MAS.xls
library(ath1121501.db); library(genefilter)

mean.mas <- read.table("Mean_Intensities_MAS.txt", header=TRUE, row.names=1, sep="\t", fill=TRUE)
sub1 <- sub("C$", "__Control", colnames(mean.mas)); sub2 <- sub("H$", "__Hypoxia", sub1) 
colnames(mean.mas) <- sub2

mean.mas <- mean.mas[rownames(mean.slr), ]
mean.mas <- as.matrix(mean.mas); rownames(mean.mas) <- mean.slr$AgI

# Keep max for duplicated genes.
mean.mas <- aggregate(mean.mas, by=list(rownames(mean.mas)), FUN=max); # Sorted automatically.
rownames(mean.mas) <- mean.mas[, 1]; mean.mas <- mean.mas[, -c(1, ncol(mean.mas)), drop=FALSE]

# Append annotations.
ann <- mean.slr$Descriptions..Canada.; names(ann) <- mean.slr$AgI
ann <- ann[!duplicated(names(ann))]
ann <- ann[rownames(mean.mas)]
mean.mas <- cbind(mean.mas, Descriptions.Canada=ann, stringsAsFactors=FALSE)
rownames(mean.mas) <- toupper(rownames(mean.mas))

# Annotation from ath1121501.db.
ann1 <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "), stringsAsFactors=FALSE)
rm1 <- grep("AFFX|s_at|a_at|x_at|r_at", rownames(ann1)); rm2 <- which(ann1[, "SYMBOL"]=="NA"); rm3 <- which(ann1$DESC=="NA")
ann2 <- ann1[-c(rm1, rm2, rm3), ]
# Remove duplicated probes and genes.
dup <- rownames(ann2)[which(duplicated(rownames(ann2)))]
ann2 <- subset(ann2, !(rownames(ann2) %in% dup))
dup1 <- ann2[, 2][which(duplicated(ann2[, 2]))]
ann2 <- subset(ann2, !(ann2[, 2] %in% dup1))

# Only keep genes with SYMBOLs.
sym <- subset(ann2, ACCNUM %in% rownames(mean.mas))
mean.mas <- mean.mas[sym$ACCNUM, ]
rownames(mean.mas) <- make.names(sym$SYMBOL)

# Log-2 transform the intensity and filter genes.
mean.mas[, -46] <- log2(mean.mas[, -46]); dim(mean.mas)                                                                                                                                                                                         
cna <- colnames(mean.mas); cna <- sub('^RT', 'rootTip', cna); cna <- sub('^R', 'root', cna); cna <- sub('^S', 'shoot', cna)                                                                                                                     
colnames(mean.mas) <- cna; df <- mean.mas                                                                                                                                                                                                       
                                                                                                                                                                                                                                                
quantile(as.matrix(df[, -46]), probs=seq(0, 1, 0.25))                                                                                                                                                                                           
ffun <- filterfun(pOverA(1/ncol(df[, -46]), 7), cv(0.40, 10000))                                                                                                                                                                                
filtered <- genefilter(df[, -46], ffun); df <- df[filtered, ]; dim(df)                                                                                                                                                                          
write.table(df, "ucr_efp_expr_ann_row_gen.txt", row.names=TRUE, col.names=TRUE, sep="\t")








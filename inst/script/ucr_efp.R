# ucr_efp_expr_ann_row_gen.txt

# Source: http://efp.ucr.edu/data/Mean_Intensities_MAS.xls
library(ath1121501.db); library(genefilter)
options(stringsAsFactors=FALSE)

mean.mas <- read.table("Mean_Intensities_MAS.txt", header=TRUE, row.names=1, sep="\t", fill=TRUE)
sub1 <- sub("C$", "__Control", colnames(mean.mas)); sub2 <- sub("H$", "__Hypoxia", sub1) 
colnames(mean.mas) <- sub2

# Annotation from ath1121501.db.
ann1 <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "), stringsAsFactors=FALSE)
rm1 <- grep("AFFX|s_at|a_at|x_at|r_at", rownames(ann1)); rm2 <- which(ann1[, "SYMBOL"]=="NA"); rm3 <- which(ann1$DESC=="NA")
ann2 <- ann1[-c(rm1, rm2, rm3), ]

pro.sym <- subset(ann2, rownames(ann2) %in% rownames(mean.mas))[, c('ACCNUM', 'SYMBOL', 'DESC')]
int <- intersect(rownames(mean.mas), rownames(pro.sym))
mean.mas <- mean.mas[int, ]; pro.sym <- pro.sym[int, ]
mean.mas <- as.matrix(mean.mas)
rownames(mean.mas) <- make.names(pro.sym$SYMBOL)

# Keep max for duplicated genes.
mean.mas <- aggregate(mean.mas, by=list(rownames(mean.mas)), FUN=max); # Sorted automatically.
rownames(mean.mas) <- mean.mas[, 1]; mean.mas <- mean.mas[, -1, drop=FALSE]

# Log-2 transform the intensity and filter genes.
mean.mas <- log2(mean.mas); dim(mean.mas)

cna <- colnames(mean.mas); cna <- sub('^RT(?!ot)', 'rootTip', cna); cna <- sub('^R', 'root', cna); cna <- sub('^S', 'shoot', cna)
colnames(mean.mas) <- cna; df <- mean.mas  

quantile(as.matrix(df), probs=seq(0, 1, 0.25))
ffun <- filterfun(pOverA(1/ncol(df), 7), cv(0.40, 10000))
filtered <- genefilter(df, ffun); df <- df[filtered, ]; dim(df)                                               

# Append annotation.
des <- pro.sym$DESC; names(des) <- make.names(pro.sym$SYMBOL)
des <- des[rownames(df)]
df <- cbind(df, annotation=des)

write.table(df, "ucr_efp_expr_ann_row_gen.txt", row.names=TRUE, col.names=TRUE, sep="\t")








\name{adj.mod}
\alias{adj.mod}
\title{
 Compute adjacency matrix and identify modules.
}
\description{
 This function is designed for large data matrix (e.g.: > 10,000 rows), which cannot be 
 processed using online version of the app. It computes the adjacency matrix and identifies 
 moudles for the input data matrix by using functions from the R package "WGCNA".
}
\usage{
adj.mod(data, type, minSize)
}
\arguments{
  \item{data}{
   The processed data matrix resulting from the function "filter.data()". }
}
  \item{type}{
   There are two options in this argument, "signed" and "unsigned". "signed" means both 
   positive and negative adjacency between genes are maintained in network module 
   identification while "unsigned" takes the absolute values of negative adjacency. Refer to
   "WGCNA" for more details.
}
  \item{minSize}{
   The minimum module size in module identification. The default value is 15. Refer to 
   "WGCNA" for more details.
}
}
\details{
 This function relies on "WGCNA". First, it computes the adjacency matrix based on the input
 data matrix, then converts the resulting adjacency matrix to a scale-free topology using a
 soft threshold (sft). If the type is "signed" sft=12, if the type is "unsigned" sft=6. 
 Next, the modules are identified using dynamic tree cutting algorithm at four sensitivity 
 levels 0, 1, 2, 3. From 0 to 3, more modules are identified but module sizes are 
 decreasing. In order to enhance the performance of this Spatial Heatmap, modules are only
 identified at level 2 and 3, since the interactive network functionality of this app works
 better on smaller modules. Refer to "WGCNA" for more details.
}
\value{
 The adjacency matrix file "adj.txt" and module assignment file "mod.txt", which are all 
 automatically saved in the directory "local_mode_result."
}
\references{
 Langfelder P and Horvath S, WGCNA: an R package for weighted correlation
 network analysis. BMC Bioinformatics 2008, 9:559 doi:10.1186/1471-2105-9-559

 Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
 Correlations and Hierarchical Clustering. Journal of Statistical Software,
 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.

 Peter Langfelder, Steve Horvath (2012). Fast R Functions for Robust
 Correlations and Hierarchical Clustering. Journal of Statistical Software,
 46(11), 1-17. URL http://www.jstatsoft.org/v46/i11/.
}
\author{
 Jianhai Zhang <jzhan067@ucr.edu; zhang.jianhai@hotmail.com>.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 data.path <- system.file("extdata/exmaple", "gene_expr_test.txt", package = "spatialHeatmap")
 exp <- filter.data(data=data.path, sep="\t", isRowGen=T, c(0, 0), cv(0.1, 10000),
 "processed_data")
 adj_mod <- adj.mod(data=exp, type="signed", minSize=20)

## The function is currently defined as
function (data, type, minSize) 
{
    library(WGCNA)
    library(flashClust)
    if (!dir.exists("local_mode_result")) 
        dir.create("local_mode_result")
    if (type == "signed") 
        sft <- 12
    if (type == "unsigned") 
        sft <- 6
    data <- t(data[, grep("__", colnames(data)), drop = F])
    adj = adjacency(data, power = sft, type = type)
    diag(adj) = 0
    tom <- TOMsimilarity(adj, TOMType = "signed")
    dissTOM = 1 - tom
    tree.hclust = flashClust(as.dist(dissTOM), method = "average")
    mcol <- NULL
    for (ds in 2:3) {
        tree <- cutreeHybrid(dendro = tree.hclust, pamStage = F, 
            minClusterSize = (minSize - 3 * ds), cutHeight = 0.99, 
            deepSplit = ds, distM = dissTOM)
        mcol <- cbind(mcol, tree$labels)
    }
    colnames(mcol) <- as.character(2:3)
    rownames(mcol) <- 1:nrow(mcol)
    write.table(adj, paste0("./local_mode_result/", "adj.txt"), sep = "\t", row.names = T, 
        col.names = T)
    write.table(mcol, paste0("./local_mode_result/", "mod.txt"), 
        sep = "\t", row.names = T, col.names = T)
    adj.mod <- list(adj = adj, mod = mcol)
    return(adj.mod)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

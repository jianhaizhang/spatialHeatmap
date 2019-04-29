data.path <- system.file("extdata/example", "root_expr_ann_row_gen.txt", package = "spatialHeatmap")

test_filter.data <- function() {
    exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), dir="./")
    checkTrue(is(exp, "data.frame"))
}

test_adj.mod <- function() {
    exp <- filter.data(data=data.path, sep="\t", isRowGen=TRUE, pOA=c(0, 0), CV=c(0.1, 10000), dir="./")
    adj_mod <- adj.mod(data=exp, type="signed", minSize=20)
    checkTrue(is(adj_mod, "list"))
}


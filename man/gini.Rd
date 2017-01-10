\name{gini}
\alias{gini}
\title{Gini Correlation Coefficient}
\usage{
gini(edges, expression, bootstrapIterations, statCutoff)
}
\description{
For each of the edges calculate the gini correlation coefficient.
If it is greater than \code{statCutoff}, also calculate a p-value by doing
\code{bootstrapIterations} random permutations of the expression data.
Returns an \R \code{data.frame} containing source ,target, gini, and pval columns
}
\examples{
edges <- read.delim("data/network.tsv")
expr <- as.matrix(read.delim("data/expression.tsv", row.names = 1))
gini(edges, expr, 10000, 0.5)
}

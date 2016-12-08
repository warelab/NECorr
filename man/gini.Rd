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
edges <- data.frame(source=c("a","a"),target=c("b","c"))
expr <- data.frame(a=c(1,2,3,4,5),b=c(1,2,3,5,6),c=c(5,4,3,2,3))
gini(edges, expr, 10000, 0.5)
}

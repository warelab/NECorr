% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorrelationGini.R
\name{cor.pair}
\alias{cor.pair}
\title{cor.pair}
\usage{
cor.pair(
  idxvec,
  GEMatrix,
  rowORcol = c("row", "col"),
  cormethod = c("GCC", "PCC", "SCC", "KCC"),
  pernum = 0,
  sigmethod = c("two.sided", "one.sided")
)
}
\arguments{
\item{idxvec}{vector of the index}

\item{GEMatrix}{Gene expression network}

\item{rowORcol}{row or columns}

\item{cormethod}{correlation type "GCC","PCC", "SCC" or "KCC"}

\item{pernum}{permutation number}

\item{sigmethod}{p-value threshold}
}
\value{
res paired correlations
}
\description{
cor.pair
}

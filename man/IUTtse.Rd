% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TissueSpecificity.R
\name{IUTtse}
\alias{IUTtse}
\title{IUTtse}
\usage{
IUTtse(data, geneTS, targetSample, expThreshold)
}
\arguments{
\item{data}{expression data}

\item{geneTS}{genes with tissue specific proprieties}

\item{targetSample}{target sample name}

\item{expThreshold}{threshold of the expression for RNAseq around 10 for microarray
could be different according to the platform}
}
\value{
y.sort
}
\description{
define the tissue excluded genes from the tissue/condition of interest
}

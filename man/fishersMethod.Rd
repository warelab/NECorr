% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{fishersMethod}
\alias{fishersMethod}
\title{fishersMethod}
\usage{
fishersMethod(x)
}
\arguments{
\item{x}{vector of p-values}
}
\value{
res
}
\description{
combine p-values Fisher method to have overall importance of
this node in the studied tissue (using code from Michael Love code)
(http://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/)
}

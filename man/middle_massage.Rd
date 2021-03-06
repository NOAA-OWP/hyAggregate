% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{middle_massage}
\alias{middle_massage}
\title{Re-index the interior of vector by threshold
Merge the interior values of a vector if they are less then the provided threshold.
Merging will look "up" and "down" the vector and merge into the smaller of the two.}
\usage{
middle_massage(x, ind, thres)
}
\arguments{
\item{x}{vector of values}

\item{ind}{current index values}

\item{thres}{threshold to evaluate x}
}
\value{
a vector of length(x) containing grouping indexes
}
\description{
Re-index the interior of vector by threshold
Merge the interior values of a vector if they are less then the provided threshold.
Merging will look "up" and "down" the vector and merge into the smaller of the two.
}

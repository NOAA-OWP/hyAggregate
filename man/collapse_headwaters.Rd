% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregation.R
\name{collapse_headwaters}
\alias{collapse_headwaters}
\title{Collapse Headwaters}
\usage{
collapse_headwaters(
  network_list,
  min_area_sqkm = 3,
  min_length_km = 1,
  verbose = TRUE,
  cache_file = NULL
)
}
\arguments{
\item{network_list}{a list containing flowpath and catchment \code{sf} objects}

\item{min_area_sqkm}{The minimum allowable size of the output hydrofabric catchments}

\item{min_length_km}{The minimum allowable length of the output hydrofabric flowlines}

\item{verbose}{should messages be emitted?}

\item{cache_file}{If not NULL results will be written to a provide path (.gpkg)}
}
\value{
a list containing flowpath and catchment \code{sf} objects
}
\description{
This function identifies small (pathlength or area) headwater catchments and
collapses them into the existing network until none remain.
Headwaters are those segments in which there are no inflows (!ID \%in\% toID).
}

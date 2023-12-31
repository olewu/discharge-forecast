% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_donors.R
\name{find_my_donors}
\alias{find_my_donors}
\title{Find donor catchments based on catchment properties}
\usage{
find_my_donors(myproperties, theirproperties, properties_to_include)
}
\arguments{
\item{myproperties}{A data table with (standardized) catchment attributes of the target catchment (typically ungauged).}

\item{theirproperties}{A data table with (standardized) catchment attributes for potential donor catchments (typically gauged).}

\item{properties_to_include}{Vector with the name of the attributes/covariates that should be included in the donor search.}
}
\value{
A data table similar to "theirproperties" with the catchments' (L2-)distance to the target catchment included. The distance is here taken across all catchments attributes mentioned in "properties_to_include".
}
\description{
Find donor catchments based on catchment properties
}
\examples{
properties_to_include=c("utm_east_z33","utm_north_z33")
myproperties=scaledprop[stat_id=="42.16.0",]
theirproperties=scaledprop[stat_id!="42.16.0",]
find_my_donors(myproperties,theirproperties,properties_to_include)
}

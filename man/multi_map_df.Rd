\name{multi_map_df}
\alias{multi_map_df}
\title{
multi_map_df
}
\description{
  mapping function (it add the possibility to map using more than one column of the data frame)
}
\usage{
multi_map_df(dfToMap, dfMap, strColsFrom, strColFromDfMap, strColToDfMap, caseSensitive=FALSE)
}

\arguments{
  \item{dfToMap}{
input data frame (that contains the columns that need to be mapped)
}
  \item{dfMap}{
data frame containing the mapping data
}
  \item{strColsFrom}{
sorted vector containing the names of the columns to be used in the input data frame for the mapping (the order of the elements in the vector defines the priority for the mapping)
}
  \item{strColFromDfMap}{
name of the column in the mapping data frame to be used as source for the mapping
}
  \item{strColToDfMap}{
name of the column in the mapping data frame to be used as target for the mapping
}
  \item{caseSensitive}{
specify whether the mapping should be case sensitive
}
}

\value{
data frame with an additional column containing the result of the mapping
}
\author{
Andrea Franceschini
}

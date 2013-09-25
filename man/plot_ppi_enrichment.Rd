\name{plot_ppi_enrichment}
\alias{\S4method{plot_ppi_enrichment}{STRINGdb}}
\title{
plot_ppi_enrichment
}
\description{
Plots a graph showing the enrichment in protein interactions in various positions of the list of genes in input.
}
\usage{
\S4method{plot_ppi_enrichment}{STRINGdb}(string_ids, file=NULL, sliceWindow = 20, edgeWindow = 140, 
                                     windowExtendedReferenceThreshold = 260, minVal=0.0000000001, title="", quiet=FALSE)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers
}
  \item{file}{
file where to save the graph as an image
}
  \item{sliceWindow}{
defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
}
  \item{edgeWindow}{
size of the window that we use to compute the enrichment (i.e. the window pvalue is computed using the proteins inside this "edgeWindow")
}
  \item{windowExtendedReferenceThreshold}{
defines the size of a window at the beginning of the list. 
                                              The enrichment will be computed always including the proteins in this window
}
  \item{title}{
title of the graph.
}
  \item{minVal}{
minimum value that the pvalue can assume in the log-scale graph. If the p-value is lower, we convert the value to this minimum value
}
  \item{quiet}{
if set to TRUE the method runs in quiet mode (turning off any output message)
}
}
\author{
Andrea Franceschini
}



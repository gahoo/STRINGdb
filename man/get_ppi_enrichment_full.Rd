\name{get_ppi_enrichment_full}
\alias{\S4method{get_ppi_enrichment_full}{STRINGdb}}
\title{
get_ppi_enrichment_full
}
\description{
Returns a vector showing the enrichment in protein interactions in various positions of the list of genes in input.
In practice, a list of 3 vectors is returned:
1) enrichment  (i.e.  enrichment computed in the window from 1 to x)
2) enrichmentWindow (i.e. enrichment computed in a sliding window of size determined by the "edgeWindow" parameters 
and the sliding steps determined by the "sliceWindow" parameter)
3) enrichmentWindowExtended  (i.e. like the enrichmentWindow, 
but it also includes an initial window of size "windowExtendedReferenceThreshold" with respect to which to compute the enrichment )
}
\usage{
\S4method{get_ppi_enrichment_full}{STRINGdb}(string_ids, sliceWindow = 20, edgeWindow  = 140, windowExtendedReferenceThreshold = 260, growingWindowLimit=NULL)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers
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
  \item{growingWindowLimit}{
threshold where to stop the computation of the enrichment
}
}
\value{
Returns a vector showing the enrichment in protein interactions in various positions of the list of genes in input.
}
\author{
Andrea Franceschini
}



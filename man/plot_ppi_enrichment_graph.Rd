\name{plot_ppi_enrichment_graph}
\alias{plot_ppi_enrichment_graph}

\title{
plot_ppi_enrichment_graph
}
\description{
Plots a graph showing the enrichment in protein interactions in various positions of the list of genes in input.
}
\usage{
plot_ppi_enrichment_graph(proteins, ppi_network, file, sliceWindow, edgeWindow, windowExtendedReferenceThreshold, minVal, title, quiet)
}

\arguments{
  \item{proteins}{
a vector of protein identifiers
}
\item{ppi_network}{
an igraph object containing the protein-protein interactions' graph.
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
\references{
Pradines JR, Farutin V, Rowley S, Dancik V. Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees. J. Comput. Biol. 2005;12:113-128.
Franceschini, A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:Nucleic Acids Res. 2013 Jan;41(Database issue)
}
\author{
Andrea Franceschini
}


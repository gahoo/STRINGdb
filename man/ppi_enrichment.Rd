\name{ppi_enrichment}
\alias{ppi_enrichment}

\title{
ppi_enrichment
}
\description{
Computes the enrichment in protein-protein interactions.
}
\usage{
ppi_enrichment(hitList, ppi_network)
}

\arguments{
  \item{hitList}{
  sorted list of proteins (from the most significant to the least significant)
}
  \item{ppi_network}{
  an igraph object containing the graph of the protein-protein interaction's network.
}
}

\value{
  \item{enrichment}{pvalue that describes the probability to get such a number of interactions by chance}
  \item{lambda}{expected number of interactions}
}
\references{
Pradines JR, Farutin V, Rowley S, Dancik V. Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees. J. Comput. Biol. 2005;12:113-128.
Franceschini, A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:Nucleic Acids Res. 2013 Jan;41(Database issue)
}
\author{
Andrea Franceschini
}


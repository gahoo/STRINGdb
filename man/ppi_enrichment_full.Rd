\name{ppi_enrichment_full}
\alias{ppi_enrichment_full}

\title{
ppi_enrichment_full
}
\description{
Compute the enrichment in protein-protein interactions of a sorted list of proteins.  
The computation is repeated at different positions in the list.  
}
\usage{
ppi_enrichment_full(hitList, ppi_network, sliceWindow, edgeWindow, windowExtendedReferenceThreshold, growingWindowLimit, quiet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hitList}{
  sorted list of proteins (from the most significant to the least significant)
}
\item{ppi_network}{
an igraph object containing the protein-protein interactions' graph.
}
\item{sliceWindow}{
defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
}
\item{edgeWindow}{
size of the window that we use to compute the enrichment (i.e. the window pvalue is computed using the proteins inside this "edgeWindow")
}
\item{windowExtendedReferenceThreshold}{
when we compute the "windowExtended" pvalue we are computing the pvalue that considers the following interactions:
   1) the interactions inside the edgeWindow (as we do with the edgeWindow pvalue)
   2) the interactions that connects the proteins in the edgeWindow with the proteins in another window at the beginning of the list (i.e. the windowExtendedReference).
windowExtendedReferenceThrehold defines the size of this windowExtendedReference window.
In this way we can compute, in a reliable way, the enrichment of a sorted list of proteins, in various positions of the list.
}
\item{growingWindowLimit}{
stop to compute the enrichment (from position 1 to position n) after growingWindowLimit proteins in the sorted list.
  (this limit speeds up the computation of the 2 other types of enrichment)
}
\item{quiet}{
if set to TRUE the method runs in quiet mode (turning off any output message)
}
}

\value{
\item{enrichment}{vector containing the enrichments in protein-protein interactions of the input list of genes (the length of the vector is length(hitList)/sliceWindow ).    }
\item{enrichmentWindow}{vector containing the enrichments in protein-protein interactions of the input list of genes (the length of the vector is length(hitList)/sliceWindow ).  The enrichment is computed considering only the proteins inside the sliding window}
\item{enrichmentWindowExtended}{vector containing the enrichments in protein-protein interactions of the input list of genes (the length of the vector is length(hitList)/sliceWindow ).
Look at the description of the windowExtendedReferenceThreshold variable
}
}
\references{
Pradines JR, Farutin V, Rowley S, Dancik V. Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees. J. Comput. Biol. 2005;12:113-128.
Franceschini, A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:Nucleic Acids Res. 2013 Jan;41(Database issue)
}
\author{
Andrea Franceschini
}


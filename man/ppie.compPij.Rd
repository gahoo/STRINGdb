\name{ppie.compPij}
\alias{ppie.compPij}

\title{
ppie.compPij
}
\description{
Compute the pvalue of protein I to interact with protein J (look at the Pradines paper referenced below)
}
\usage{
ppie.compPij(degI, degJ, edgeNum)
}

\arguments{
  \item{degI}{
degree of protein I
}
 \item{degJ}{
degree of protein J
}
 \item{edgeNum}{
total number of edges of the entire graph
}
}
\value{
return the pvalue of protein I to interact with protein J (look at the Pradines paper referenced below)
}
\references{
Pradines JR, Farutin V, Rowley S, Dancik V. Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees. J. Comput. Biol. 2005;12:113-128.
Franceschini, A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:Nucleic Acids Res. 2013 Jan;41(Database issue)
}
\author{
Andrea Franceschini
}

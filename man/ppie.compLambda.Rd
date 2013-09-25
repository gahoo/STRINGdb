\name{ppie.compLambda}
\alias{ppie.compLambda}

\title{
ppie.compLambda
}
\description{
compute the number of expected interactions between a set of proteins.
}
\usage{
ppie.compLambda(degrees, edgeNum)
}

\arguments{
  \item{degrees}{
vector containing the degrees of the nodes of a set of proteins
}
  \item{edgeNum}{
total number of edges of the entire graph
}
}

\value{
number of expected interactions
}
\references{
Pradines JR, Farutin V, Rowley S, Dancik V. Analyzing protein lists with large networks: edge-count probabilities in random graphs with given expected degrees. J. Comput. Biol. 2005;12:113-128.
Franceschini, A et al. (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:Nucleic Acids Res. 2013 Jan;41(Database issue)
}
\author{
Andrea Franceschini
}
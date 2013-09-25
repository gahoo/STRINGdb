\name{mp}
\alias{\S4method{mp}{STRINGdb}}
\title{
mp
}
\description{
Maps the gene identifiers of the input vector to STRING identifiers (using a take first approach).
It returns a vector with the STRING identifiers of the mapped proteins.
}
\usage{
\S4method{mp}{STRINGdb}(protein_aliases)
}
\arguments{
  \item{protein_aliases}{
vector of protein aliases that we want to convert to STRING identifiers
}

}
\value{
It returns a vector with the STRING identifiers of the mapped proteins.
}
\author{
Andrea Franceschini
}



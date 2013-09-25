\name{get_ppi_enrichment}
\alias{\S4method{get_ppi_enrichment}{STRINGdb}}
\title{
get_ppi_enrichment
}
\description{
Returns a pvalue representing the enrichment in interactions of the list of proteins 
          (i.e. the probability to obtain such a number of interactions by chance).
}
\usage{
\S4method{get_ppi_enrichment}{STRINGdb}(string_ids)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers
}
}
\value{
Returns a pvalue representing the enrichment in interactions of the list of proteins  (i.e. the probability to obtain such a number of interactions by chance).
}
\author{
Andrea Franceschini
}



\name{get_term_proteins}
\alias{\S4method{get_term_proteins}{STRINGdb}}
\title{
get_term_proteins
}
\description{
Returns the proteins annotated to belong to a given term.
}
\usage{
\S4method{get_term_proteins}{STRINGdb}(term_ids, string_ids=NULL, enableIEA=TRUE)
}
\arguments{
  \item{term_ids}{
vector of terms
}
  \item{string_ids}{
a vector of STRING identifiers. If the variable is set, the method returns only the proteins that are present in this vector.
}
  \item{enableIEA}{
whether to consider also Electronic Inferred Annotations
}
}
\value{
Returns the proteins annotated to belong to a given term.
}
\author{
Andrea Franceschini
}



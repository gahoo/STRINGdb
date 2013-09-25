\name{get_enrichment}
\alias{\S4method{get_enrichment}{STRINGdb}}
\title{
get_enrichment
}
\description{
Returns the enrichment in pathways of the vector of STRING proteins that is given in input.
}
\usage{
\S4method{get_enrichment}{STRINGdb}(string_ids, category = "Process", methodMT = "fdr", iea = TRUE)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers.
}
  \item{category}{
category for which to compute the enrichment (i.e. "Process", "Component", "Function", "KEGG", "Pfam", "InterPro"). The default category is "Process".
}
  \item{methodMT}{
method to be used for the multiple testing correction. (i.e.  "fdr", "bonferroni"). The default is "fdr".
}
  \item{iea}{
specify whether you also want to use electronic inference annotations
}
}
\value{
Data frame containing the enrichment in pathways of the vector of STRING proteins that is given in input.
}
\author{
Andrea Franceschini
}



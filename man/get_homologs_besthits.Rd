\name{get_homologs_besthits}
\alias{\S4method{get_homologs_besthits}{STRINGdb}}
\title{
get_homologs_besthits
}
\description{
Returns the best blast hits x species of the given input identifiers.
}
\usage{
\S4method{get_homologs_besthits}{STRINGdb}(string_ids, symbets = FALSE, target_species_id = NULL, bitscore_threshold=NULL)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers.
}
  \item{target_species_id}{
 NCBI taxonomy identifier of the species to query for homologs (the species must be present in the STRING database)
}
  \item{bitscore_threshold}{
threshold on the bitscore of the blast alignment.
}
  \item{symbets}{
specify whether you want only symmetrical best hits
}
}
\value{
Data frame containing the best blast hits x species of the given input identifiers.
}
\author{
Andrea Franceschini
}



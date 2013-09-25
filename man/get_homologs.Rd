\name{get_homologs}
\alias{\S4method{get_homologs}{STRINGdb}}
\title{
get_homologs
}
\description{
Returns the homologs of the given input identifiers that are present in the given target_species_id.
}
\usage{
\S4method{get_homologs}{STRINGdb}(string_ids, target_species_id, bitscore_threshold=NULL)
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
}
\value{
Data frame containing the homologs of the given input identifiers and that are present in the given target_species_id.
}
\author{
Andrea Franceschini
}



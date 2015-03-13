\name{remove_homologous_interactions}
\alias{\S4method{remove_homologous_interactions}{STRINGdb}}
\title{
remove_homologous_interactions
}
\description{
With this method it is possible to remove the interactions that are composed by a pair of homologous/similar proteins, having a similarity bitscore between each other higher than a threshold.
}
\usage{
\S4method{remove_homologous_interactions}{STRINGdb}(interactions_dataframe, bitscore_threshold = 60)
}
\arguments{
  \item{interactions_dataframe}{
 a data frame contaning the sorted interactions to be benchmarked. The data frame should have the following column names: proteinA, proteinB, score
}
  \item{bitscore_threshold}{
filter out pairs of homologous proteins, having a similarity bitscore higher than this parameter
}
}
\value{
interactions data frame where the homologous pairs have been removed, from the input interactions' data frame
}
\author{
Andrea Franceschini
}

\name{benchmark_ppi}
\alias{\S4method{benchmark_ppi}{STRINGdb}}
\title{
benchmark_ppi
}
\description{
benchmark a list of protein-protein interactions using pathways (e.g. KEGG). 
The function outputs a table where the interactions are mapped to KEGG and the number of TPs and FPs are counted.
}
\usage{
\S4method{benchmark_ppi}{STRINGdb}(interactions_dataframe, pathwayType = "KEGG", max_homology_bitscore = 60, precision_window=400, exclude_pathways="blacklist")
}
\arguments{
  \item{interactions_dataframe}{
 a data frame contaning the sorted interactions to be benchmarked. The data frame should have the following column names: proteinA, proteinB, score
}
  \item{pathwayType}{
category to use to benchmark the interactions (default KEGG)
}
  \item{max_homology_bitscore}{
filter out pairs of homologous proteins, having a similarity bitscore higher than this parameter
}
  \item{precision_window}{
define the size of the precision window (i.e. the window used to scan the sorted interactions data frame in order to compute the benchmark precision).
At the beginning and at the end of the sorted interactions data frame, the window is automatically enlarged (at the beginning) and reduced (at the end)
}
  \item{exclude_pathways}{
Exlude the terms that should not be used for benchmarking. If this parameter is set to "blacklist", a black list to be used mantained by our group is automatically downloaded from our servers (otherwise it is possible to specify a vector with the terms that have to be excluded).
}
}
\value{
interactions data frame where the interactions are mapped to KEGG and the number of TPs and FPs are counted.
}
\author{
Andrea Franceschini
}

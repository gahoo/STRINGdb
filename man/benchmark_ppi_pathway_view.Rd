\name{benchmark_ppi_pathway_view}
\alias{\S4method{benchmark_ppi_pathway_view}{STRINGdb}}
\title{
benchmark_ppi_pathway_view
}
\description{
Takes in input the results of the benchmark_ppi function, and constructs a new table that provides a view at the pathway level 
    (i.e. it list all the pathways to which the interactions belong)
}
\usage{
\S4method{benchmark_ppi_pathway_view}{STRINGdb}(benchmark_ppi_data_frame, precision_threshold=0.2, pathwayType = "KEGG")
}
\arguments{
  \item{benchmark_ppi_data_frame}{
data frame that comes out from the benchmark_ppi function.
}
  \item{precision_threshold}{
threshold that specify where to stop taking considering the interactions in the sorted input data frame (the list is scanned until the precision goes below this value).
}
  \item{pathwayType}{
the pathway category to use (KEGG by default)
}
}
\value{
data frame containing the pathways (i.e. terms) of the input proteins. 
Several parameters are reported: 
pathway coverage (i.e. number of interactions in the list belonging to the pathway / maximum number of interactions between the pathay proteins (i.e. proteins * (proteins-1)/2 ) )
total_representation (i.e. number of interactions in the list belonging to the pathway / size of the input data frame )
}
\author{
Andrea Franceschini
}

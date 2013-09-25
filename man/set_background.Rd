\name{set_background}
\alias{\S4method{set_background}{STRINGdb}}
\title{
set_background
}
\description{
With this method you can specify a vector of proteins to be used as background. 
The network is reloaded and only the proteins that are present in the background vector are inserted in the graph.  
Besides, the background is taken in consideration for all the enrichment statistics.
}
\usage{
\S4method{set_background}{STRINGdb}(background_vector )
}
\arguments{
  \item{background_vector}{
vector of STRING protein identifiers
}
}
\author{
Andrea Franceschini
}



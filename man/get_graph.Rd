\name{get_graph}
\alias{\S4method{get_graph}{STRINGdb}}
\title{
get_graph
}
\description{
Return an igraph object with the STRING network (for information about iGraph visit http://igraph.sourceforge.net)
}
\references{
  Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. 
  http://igraph.sf.net
}
\seealso{
In order to simplify the most common tasks, we do also provide convenient functions that wrap some iGraph functions.
  get_interactions(string_ids)   # returns the interactions in between the input proteins
  get_neighbors(string_ids)      # Get the neighborhoods of a protein (or of a vector of proteins) that is given in input.
  get_subnetwork(string_ids)     # returns a subgraph from the given input proteins
}
\usage{
\S4method{get_graph}{STRINGdb}()
}
\value{
igraph object with the STRING network
}
\author{
Andrea Franceschini
}



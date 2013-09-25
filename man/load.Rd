\name{load}
\alias{\S4method{load}{STRINGdb}}
\title{
load
}
\description{
Downloads and returns the STRING network (the network is set also in the graph variable of the STRING_db object).

  
It makes use of the variables:
    "backgroundV"         vector containing STRING identifiers to be used as background 
                            (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
    "score_threshold"     STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)
}
\usage{
\S4method{load}{STRINGdb}()
}
\value{
STRING network (i.e. an iGraph object.  For info look to http://igraph.sourceforge.net)
}
\author{
Andrea Franceschini
}



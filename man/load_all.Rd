\name{load_all}
\alias{\S4method{load_all}{STRINGdb}}
\title{
load_all
}
\description{
Force download and loading of all the files (so that you can later store the object on the hard disk if you like).
It makes use of the variables:
"backgroundV"       vector containing STRING identifiers to be used as background 
                    (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
"score_threshold"   STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)
}
\usage{
\S4method{load_all}{STRINGdb}()
}
\author{
Andrea Franceschini
}



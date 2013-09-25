\name{get_png}
\alias{\S4method{get_png}{STRINGdb}}
\title{
get_png
}
\description{
Returns a png image of a STRING protein network with the given identifiers.
}
\usage{
\S4method{get_png}{STRINGdb}(string_ids, required_score=NULL, network_flavor="evidence", file=NULL, payload_id=NULL)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers.
}
  \item{required_score}{
minimum STRING combined score of the interactions 
(if left NULL we get the combined score of the object, which is 400 by default).
}
  \item{network_flavor}{
specify the flavor of the network ("evidence", "confidence" or "actions".  default "evidence").
}
  \item{file}{
file where to save the image
}
  \item{payload_id}{
identifier of the payload
}
}
\value{
Returns a png image of a STRING protein network with the given identifiers.
}
\author{
Andrea Franceschini
}



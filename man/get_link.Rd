\name{get_link}
\alias{\S4method{get_link}{STRINGdb}}
\title{
get_link
}
\description{
Returns a short link to the network page of our STRING website that shows the protein interactions between the given identifiers.
}
\usage{
\S4method{get_link}{STRINGdb}(string_ids, required_score=NULL, network_flavor="evidence", payload_id = NULL)
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
  \item{payload_id}{
an identifier of payload data on the STRING server (see method post_payload for additional informations)
}
}
\value{
short link to the network page of our STRING website that shows the protein interactions between the input identifiers.
}
\author{
Andrea Franceschini
}



\name{plot_network}
\alias{\S4method{plot_network}{STRINGdb}}
\title{
plot_network
}
\description{
Plots an image of the STRING network with the given proteins.
}
\usage{
\S4method{plot_network}{STRINGdb}(string_ids, payload_id=NULL, required_score=NULL, add_link=TRUE, add_summary=TRUE)
}
\arguments{
  \item{string_ids}{
a vector of STRING identifiers
}
  \item{payload_id}{
an identifier of payload data on the STRING server (see method post_payload for additional informations)
}
  \item{required_score}{
a threshold on the score that overrides the default score_threshold, that we use only for the picture
}
  \item{add_link}{
parameter to specify whether you want to generate and add a short link to the relative page in STRING. 
As default this option is active but we suggest to deactivate it in case one is generating many images (e.g. in a loop). 
Deactivating this option avoids to generate and store a lot of short-urls on our server.
}
  \item{add_summary}{
parameter to specify whether you want to add a summary text to the picture. This summary includes a p-value and the number of proteins/interactions.
}
}
\author{
Andrea Franceschini
}



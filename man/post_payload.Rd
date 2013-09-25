\name{post_payload}
\alias{\S4method{post_payload}{STRINGdb}}
\title{
post_payload
}
\description{
Posts the input to STRING and returns an identifier that you can use to access the payload when you enter in our website.
}
\usage{
\S4method{post_payload}{STRINGdb}(stringIds, colors=NULL, comments=NULL, links=NULL, iframe_urls=NULL, logo_imgF=NULL, legend_imgF=NULL )
}
\arguments{
  \item{stringIds}{
vector of STRING identifiers.
}
  \item{colors}{
vector containing the colors to use for a every STRING identifier ( the order of the elements must match those in the string_ids vector)
}
  \item{comments}{
vector containing the comments to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
}
  \item{links}{
vector containing the links to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
}
  \item{iframe_urls}{
vector containing the urls of the iframes to use for every STRING identifier 
( the order of the elements must match those in the string_ids vector).
}
  \item{logo_imgF}{
path to a file containing the logo image to be display in the STRING website
}
  \item{legend_imgF}{
path to a file containing a legend image to be display in the STRING website
}
}
\value{
identifier of the payload.
}
\author{
Andrea Franceschini
}



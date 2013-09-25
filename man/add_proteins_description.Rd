\name{add_proteins_description}
\alias{\S4method{add_proteins_description}{STRINGdb}}
\title{
add_proteins_description
}
\description{
Add description coluns to the proteins that are present in the data frame given in input. 
The data frame must contain a column named "STRING_id".
}
\usage{
\S4method{add_proteins_description}{STRINGdb}(screen)
}
\arguments{
  \item{screen}{
Dataframe containing the results of the experiment (e.g. the analyzed results of a microarray or RNAseq experiment)
}
}
\value{
returns the same dataframe given in input with an additional columns containing a description of the proteins.
}
\author{
Andrea Franceschini
}



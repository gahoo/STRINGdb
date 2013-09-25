\name{map}
\alias{\S4method{map}{STRINGdb}}
\title{
map
}
\description{
Maps the gene identifiers of the input dataframe to STRING identifiers.
It returns the input dataframe with the "STRING_id" additional column.
}
\usage{
\S4method{map}{STRINGdb}(my_data_frame, my_data_frame_id_col_names, takeFirst=TRUE, removeUnmappedRows=FALSE, quiet=FALSE)
}
\arguments{
  \item{my_data_frame}{
data frame provided as input.
}
  \item{my_data_frame_id_col_names}{
vector contatining the names of the columns of "my_data_frame" that have to be used for the mapping.
}
  \item{takeFirst}{
 boolean indicating what to do in case of multiple STRING proteins that map to the same name. 
 If TRUE, only the first of those is taken. Otherwise all of them are used. (default TRUE)
}
  \item{removeUnmappedRows}{
remove the rows that cannot be mapped to STRING 
(by default those lines are left and their STRING_id is set to NA).
}
  \item{quiet}{
Setting this variable to TRUE we can avoid printing the warning relative to the unmapped values.
}

}
\value{
Returns the dataframe that is given in input with the "STRING_id" additional column.
}
\author{
Andrea Franceschini
}



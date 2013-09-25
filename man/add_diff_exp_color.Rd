\name{add_diff_exp_color}
\alias{\S4method{add_diff_exp_color}{STRINGdb}}
\title{
add_diff_exp_color
}
\description{
Take in input a dataframe containing a logFC column that reports the logarithm of the difference in expression level.
Add a "color" column to the data frame such that strongly downregulated genes are colored in green and strong upregulated genes are in red.
When the down or up-regulation is instead weak the intensity of the color gets weaker as well, accordingly.
}
\usage{
\S4method{add_diff_exp_color}{STRINGdb}(screen, logFcColStr="logFC" )
}
\arguments{
  \item{screen}{
Dataframe containing the results of the experiment (e.g. the analyzed results of a microarray or RNAseq experiment)
}
  \item{logFcColStr}{
name of the colum that contains the logFC of the expression 
}
}
\value{
vector containing the colors
}
\author{
Andrea Franceschini
}



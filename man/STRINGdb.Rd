\name{STRINGdb}
\alias{STRINGdb}

\title{
STRINGdb (an R interface to http://string-db.org)
}
\description{
The R package STRINGdb provides a convenient interface to the STRING protein-protein interactions database for the R/bioconductor users.
Please look at the manual/vignette to get additional informationd and examples on how to use the package.
STRING is a database of known and predicted protein-protein interactions. It contains information from numerous 
sources, including experimental repositories, computational prediction methods and public text collections. Each interaction is associated with a combined confidence score that integrates the various evidences. 
STRING is regularly updated , the latest version 9.05 contains information on 5 millions proteins from more than 1100 species.
The STRING web interface is freely accessible at: http://string-db.org/
}
\usage{STRINGdb(...)}
\arguments{
\item{\dots}{
description
}
}

\source{
http://string-db.org
}
\references{
Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29'.
}
\examples{
library(STRINGdb)
data(diff_exp_example1)

# create a new STRING_db object
string_db <- STRINGdb$new()

# map to STRING
example1_mapped = string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )

# get the best 200 hits
hits = example1_mapped$STRING_id[1:200]  

# plot the STRING network png 
string_db$plot_network( hits )

# plot a protein-protein enrichment graph of the best 1000 hits in order to see how the ppi signal is distributed along the sorted list
string_db$plot_ppi_enrichment( example1_mapped$STRING_id[1:1000] )


######## use the "halo" mechanism in STRING to distinguish between down and up-regulated genes  ########

# filter by p-value and add a color column (i.e. green down-regulated gened and red for up-regulated genes)
example1_mapped_pval05 = string_db$add_diff_exp_color( subset(example1_mapped, pvalue<0.05), logFcColStr="logFC" )    

# post payload information to the STRING server
payload_id = string_db$post_payload( example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color )

# display a STRING network png with the "halo"
string_db$plot_network( hits, payload_id=payload_id )


######## compute enrichment in GO annotations  ########
enrichmentGO = string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
enrichmentKEGG = string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
head(enrichmentGO, n=7)
head(enrichmentKEGG, n=7)



########  Use STRING in order to get homologous proteins in other species ########

# get the reciprocal best hits of the following protein in all the STRING species
string_db$get_homologs_besthits(c("9606.ENSP00000365757"), symbets = TRUE)

# get the homologs of the following two proteins in the mouse (i.e. species_id=10090)
string_db$get_homologs(c("9606.ENSP00000365757", "9606.ENSP00000352336"), target_species_id=10090, bitscore_threshold=60)



########  Retrieve protein-protein interactions ########

# get the neighbors (in the ppi graph) of a given protein(s)
string_db$get_neighbors(c("9606.ENSP00000365757", "9606.ENSP00000352336"))

# get the interactions in between the input proteins
string_db$get_interactions(c("9606.ENSP00000365757", "9606.ENSP00000352336"))



######## Find clusters of interactions ########

# get clusters
clustersList = string_db$get_clusters(example1_mapped$STRING_id[1:600])

# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}


}
\keyword{datasets}

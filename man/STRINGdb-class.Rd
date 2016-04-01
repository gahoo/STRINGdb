\docType{class}
\name{STRINGdb-class}
\Rdversion{1.1}
\alias{STRINGdb-class}
\alias{STRINGdb}

\title{Class \code{"STRINGdb"}}
\description{
The R package STRINGdb provides a convenient interface to the STRING protein-protein interactions database for the R/bioconductor users.
Please look at the manual/vignette to get additional informationd and examples on how to use the package.
STRING is a database of known and predicted protein-protein interactions. It contains information from numerous 
sources, including experimental repositories, computational prediction methods and public text collections. Each interaction is associated with a combined confidence score that integrates the various evidences. 
STRING is regularly updated , the latest version 9.05 contains information on 5 millions proteins from more than 1100 species.
The STRING web interface is freely accessible at: http://string-db.org/
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
Franceschini, A (2013). STRING v9.1: protein-protein interaction networks, with increased coverage and integration. In:'Nucleic Acids Res. 2013 Jan;41(Database issue):D808-15. doi: 10.1093/nar/gks1094. Epub 2012 Nov 29'.
}
\author{
Andrea Franceschini
}


\seealso{
http://stitch-db.org
}
\examples{
showClass("STRINGdb")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{annotations}:}{Object of class \code{data.frame} ~~ }
    \item{\code{annotations_description}:}{Object of class \code{data.frame} ~~ }
    \item{\code{graph}:}{Object of class \code{igraph} ~~ }
    \item{\code{proteins}:}{Object of class \code{data.frame} ~~ }
    \item{\code{speciesList}:}{Object of class \code{data.frame} ~~ }
    \item{\code{species}:}{Object of class \code{numeric} ~~ }
    \item{\code{version}:}{Object of class \code{character} ~~ }
    \item{\code{input_directory}:}{Object of class \code{character} ~~ }
    \item{\code{backgroundV}:}{Object of class \code{vector} ~~ }
    \item{\code{score_threshold}:}{Object of class \code{numeric} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{\code{set_background(background_vector)}:}{ ~~ }
    \item{\code{post_payload(stringIds, colors, comments, links, iframe_urls, logo_imgF, legend_imgF)}:}{ ~~ }
    \item{\code{plot_network(string_ids, payload_id, required_score)}:}{ ~~ }
    \item{\code{plot_ppi_enrichment(string_ids, file, sliceWindow, edgeWindow, windowExtendedReferenceThreshold, minVal, title)}:}{ ~~ }
    \item{\code{map(my_data_frame, my_data_frame_id_col_names, takeFirst, removeUnmappedRows, quiet)}:}{ ~~ }
    \item{\code{load()}:}{ ~~ }
    \item{\code{get_term_proteins(term_ids, string_ids, enableIEA)}:}{ ~~ }
    \item{\code{get_summary(string_ids)}:}{ ~~ }
    \item{\code{get_subnetwork(string_ids)}:}{ ~~ }
    \item{\code{get_ppi_enrichment_full(string_ids, sliceWindow, edgeWindow, windowExtendedReferenceThreshold, growingWindowLimit)}:}{ ~~ }
    \item{\code{get_ppi_enrichment(string_ids)}:}{ ~~ }
    \item{\code{get_proteins()}:}{ ~~ }
    \item{\code{get_png(string_ids, required_score, network_flavor, file, payload_id)}:}{ ~~ }
    \item{\code{get_neighbors(string_ids)}:}{ ~~ }
    \item{\code{get_link(string_ids, required_score, network_flavor, payload_id)}:}{ ~~ }
    \item{\code{get_interactions(string_ids)}:}{ ~~ }
    \item{\code{get_homologs_besthits(string_ids, symbets, target_species_id, bitscore_threshold)}:}{ ~~ }
    \item{\code{get_homologs(string_ids, target_species_id, bitscore_threshold)}:}{ ~~ }
    \item{\code{get_graph()}:}{ ~~ }
    \item{\code{get_enrichment(string_ids, category, methodMT, iea)}:}{ ~~ }
    \item{\code{get_clusters(string_ids, algorithm)}:}{ ~~ }
    \item{\code{get_annotations_desc()}:}{ ~~ }
    \item{\code{get_annotations()}:}{ ~~ }
    \item{\code{load_all()}:}{ ~~ }
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{add_proteins_description(screen)}:}{ ~~ }
    \item{\code{add_diff_exp_color(screen, logFcColStr)}:}{ ~~ }
    \item{\code{show()}:}{ ~~ }
  }
}

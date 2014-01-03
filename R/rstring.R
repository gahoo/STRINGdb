
setOldClass("igraph")
STRINGdb <- setRefClass("STRINGdb",
    fields = list(
      annotations="data.frame",
      annotations_description="data.frame",
      graph="igraph",
      proteins="data.frame",
      aliases_tf="data.frame",
      #speciesList="data.frame",
      species="numeric",
      version="character",    
      input_directory="character",
      backgroundV = "vector",
      score_threshold = "numeric"
    ),
    methods = list(
      
      
      
      add_diff_exp_color = function(screen, logFcColStr="logFC" ){
'
Description:
  Take in input a dataframe containing a logFC column that reports the logarithm of the difference in expression level.
  Add a "color" column to the data frame such that strongly downregulated genes are colored in green and strong upregulated genes are in red.
  When the down or up-regulation is instead weak the intensity of the color gets weaker as well, accordingly.

Author(s):
   Andrea Franceschini
'
        screen_pval05_pos = subset(screen, as.matrix(screen[logFcColStr])[, 1] > 0)
        k=exp(screen_pval05_pos[,logFcColStr])/exp(max(screen_pval05_pos[,logFcColStr]))
        screen_pval05_pos_col = data.frame(screen_pval05_pos, color = color.scale(k,c(1,1),c(1,0),c(1,0)) )
        screen_pval05_neg = subset(screen, as.matrix(screen[logFcColStr])[, 1] <= 0)
        k=exp(-screen_pval05_neg[,logFcColStr])/exp(-min(screen_pval05_neg[,logFcColStr]))
        screen_pval05_neg_col = data.frame(screen_pval05_neg, color = color.scale(k,c(1,0),c(1,1),c(1,0)) )
        return(rbind(screen_pval05_pos_col, screen_pval05_neg_col))
      },
      
      
      
      
      add_proteins_description = function(screen){
        '
Description:
  Add description coluns to the proteins that are present in the data frame given in input. 
  The data frame must contain a column named "STRING_id".

Input parameters:
      "screen"  a data frame having a "STRING_id" column 

Author(s):
   Andrea Franceschini
'
        
        proteinsDf2 = get_proteins()
        proteinsDf3 = merge(screen, proteinsDf2, by.x="STRING_id", by.y="protein_external_id", all.x=TRUE, sort=FALSE)
        return(proteinsDf3)
      },
      
      
      
      initialize=function(...) {
        callSuper(...)
        if(length(species)==0) {
          cat("WARNING: You didn't specify a species. Hence we will set 9606 (Homo Sapiens) as your species.\n")
          species <<- 9606                        
        }
        if(length(version)==0) {
          cat("WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING 9_05.\n")
          version <<- "9_05"                        
          current_version = as.character(read.table(url("http://string.uzh.ch/permanent/string/current_version"))$V1)
          if(current_version != version) cat("WARNING: A new version of the STRING R plugin has been released: we suggest you to install the latest version from Bioconductor.\n")
        }else{
          valid_versions = read.table(url("http://string.uzh.ch/permanent/string/r_lib_valid_versions"))$V1
          if(! (version %in% valid_versions)) {
            cat("ERROR: The version ", version, "is not valid. \nPlease use one of the following versions:\n")
            cat(as.character(valid_versions), sep=", ")
            stop()
          }
        }
        if(input_directory=="" || is.null(input_directory) || length(input_directory)==0) input_directory<<-tempdir()
        if(input_directory=="" || is.null(input_directory) || length(score_threshold)==0 || score_threshold<1) score_threshold <<- 1
      },
                          
      
      load_all = function(){
'
Description:
  Force download and loading of all the files (so that you can later store the object on the hard disk if you like)
  It makes use of the variables:
  "backgroundV"       vector containing STRING identifiers to be used as background 
                      (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
  "score_threshold"   STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)

Author(s):
   Andrea Franceschini
'
        
        x1 = get_annotations()
        x2 = get_annotations_desc()  
        x3 = get_proteins()
        #x4 = get_species()
        get_aliases() 
        load()
      },
      
      
      get_aliases = function(){
'
Description:
  Loads and returns STRING alias table. 

Author(s):
   Andrea Franceschini
'        
        temp = downloadAbsentFileSTRING(paste("http://string.uzh.ch/permanent/string/", version, "/protein_aliases/", species, "__protein_aliases.tsv.gz", sep=""), oD=input_directory)
        downloadAbsentFileSTRING(paste("http://string.uzh.ch/permanent/string/", version, "/protein_aliases/", species, "__protein_aliases_tf.tsv.gz", sep=""), oD=input_directory)
        aliasDf <- read.table(temp, sep = "\t", header=TRUE, quote="", stringsAsFactors=FALSE, fill = TRUE)
        aliasDf = renameColDf(aliasDf, "protein_id", "STRING_id")
        return(aliasDf)
      },
      
      
      get_annotations = function(){
'
Description:
  Loads and returns STRING annotations (i.e. GO annotations, KEGG pathways, domain databases). 
  The annotations are stored in the "annotations" variable.

Author(s):
   Andrea Franceschini
'
        
        if(nrow(annotations)==0){
          temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/enrichment_annotations/annotations/annotations_", species, ".tsv.gz", sep=""), oD=input_directory)
          ann = read.table(temp, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=FALSE)
          ann = renameColDf(ann, "V1", "STRING_id")
          ann = renameColDf(ann, "V2", "term_id")
          ann = renameColDf(ann, "V3", "category")
          annotations <<- renameColDf(ann, "V4", "type")
        }
        return(annotations)
      },
      
      
      get_annotations_desc = function(){
'
Description:
  Returns a data frame with the description of every STRING annotation term 
    (it downloads and caches the information the first time that is called).

Author(s):
   Andrea Franceschini
'
        
        if(nrow(annotations_description)==0){
          temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/enrichment_annotations/annotations_description/annotations_description.tsv.gz", sep=""), oD=input_directory)
          annotations_description <<- read.table(temp, sep="\t", header=TRUE, quote="", fill=TRUE, stringsAsFactors=FALSE)
        }
        return(annotations_description)
      },
      
      
      get_bioc_graph = function(){
'
Description:
  Returns the interaction graph as an object of the graph package in Bioconductor

Author(s):
   Andrea Franceschini
'
        print(".... please wait about 10 minutes to create the graph ....")
        if("graph" %in% installed.packages()){
          library(graph)
          if(is.null(graph)) load()
          temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/protein_links/",
                                          species, "__protein_links.tsv.gz", sep=""), oD=input_directory)
          PPI <- read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
          
          PPIselected = PPI
          if(length(score_threshold)!=0) PPIselected <- PPI[PPI$combined_score >= score_threshold,]
          if(!is.null(backgroundV)){
            PPIselected <- PPIselected[PPIselected$protein1 %in% backgroundV,]
            PPIselected <- PPIselected[PPIselected$protein2 %in% backgroundV,]
          }
          
          nel_graph <- ftM2graphNEL(cbind(PPIselected$protein1, PPIselected$protein2), W=PPIselected$combined_score)
          nel_graph2 <- ugraph(nel_graph)
          return(nel_graph2)
          #return(igraph.to.graphNEL(graph))
        }else{
          cat("ERROR: In order to run this function you must install the \"graph\" package from CRAN. \n Please install the package and run this function again.\n")
        }
      },
      
      
      get_clusters = function(string_ids, algorithm="fastgreedy"){
'
Description:
  Returns a list of clusters of interacting proteins.

References:
  See the iGraph (http://igraph.sourceforge.net/) documentation for additional information on the algorithms.
  Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. 
  http://igraph.sf.net

Input parameters:
    "string_ids"    a vector of STRING identifiers.
    "algorithm"     algorithm to use for the clustering (fastgreedy, walktrap, spinglass, edge.betweenness)

Author(s):
   Andrea Franceschini
'
        
        if(is.null(graph)) load()
        if(algorithm=="fastgreedy") fgreedy<-fastgreedy.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)
        if(algorithm=="walktrap") fgreedy<-walktrap.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)
        if(algorithm=="spinglass") fgreedy<-spinglass.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)
        if(algorithm=="edge.betweenness") fgreedy<-edge.betweenness.community(get_subnetwork(string_ids),merges=TRUE, modularity=TRUE)

        memb = membership(fgreedy)
        clusters = NULL
        for(i in 1:max(memb)){
          clusters[[i]] = names(membership(fgreedy)[membership(fgreedy)==i])
        }
        return(clusters)
      },
      
      
      
      get_enrichment = function(string_ids, category = "Process", methodMT = "fdr", iea = TRUE){
'
Description:
  Returns the enrichment in pathways of the vector of STRING proteins that is given in input.

Input parameters:
  "string_ids"    a vector of STRING identifiers.
  "category"      category for which to compute the enrichment (i.e. "Process", "Component", "Function", "KEGG", "Pfam", "InterPro".  default "Process")
  "methodMT"      method to be used for the multiple testing correction. (i.e.  "fdr", "bonferroni"..  default "fdr")

Author(s):
   Andrea Franceschini
'
        ann = get_annotations()
        temp_category = category
        ann = subset(ann, category==temp_category)
        if(!is.null(backgroundV)) ann = subset(ann, STRING_id %in% backgroundV )
        if(!iea) ann = subset(ann, type != "IEA")
        dfcount = suppressWarnings(sqldf('select term_id, count(STRING_id) as proteins from ann group by term_id', stringsAsFactors=FALSE))
        annHits = subset(ann, STRING_id %in% string_ids)
        dfcountHits = suppressWarnings(sqldf('select term_id, count(STRING_id) as hits from annHits group by term_id', stringsAsFactors=FALSE))
        dfcountMerged = merge(dfcount, dfcountHits, by.x="term_id", by.y="term_id", all.x=TRUE)
        dfcountMerged2 = data.frame(dfcountMerged, n = length(unique(ann$STRING_id)) - dfcountMerged$proteins, k = length(unique(annHits$STRING_id)))
        dfcountMerged3 = data.frame(dfcountMerged2, pvalue= phyper(dfcountMerged2$hits-1, dfcountMerged2$proteins, dfcountMerged2$n, dfcountMerged2$k, FALSE))
        dfcountMerged4 = dfcountMerged3
        if(!is.null(methodMT)) dfcountMerged4 = data.frame(dfcountMerged3, pvalue_fdr = p.adjust(dfcountMerged3$pvalue, method=methodMT, n=nrow(subset(dfcountMerged3, !is.na(pvalue))) ))
        dfcountMerged4 = subset(dfcountMerged4, !is.na(pvalue))
        dfcountMerged4 = delColDf(dfcountMerged4, "n")
        dfcountMerged4 = delColDf(dfcountMerged4, "k")
        annDesc = get_annotations_desc()
        dfcountMerged5 = arrange(merge(dfcountMerged4, annDesc, by.x="term_id", by.y="term_id", all.x=TRUE), pvalue)
        return(dfcountMerged5)
      },
      
      
      get_graph = function(){
'
Description:
  Return an igraph object with the entire STRING network. 
  We invite the user to use the functions of the iGraph package to conveniently 
  search/analyze the network.
  
References:
  Csardi G, Nepusz T: The igraph software package for complex network research, 
  InterJournal, Complex Systems 1695. 2006. 
  http://igraph.sf.net

See Also:
  In order to simplify the most common tasks, we do also provide convenient functions 
  that wrap some iGraph functions.
  get_interactions(string_ids)   # returns the interactions in between the input proteins
  get_neighbors(string_ids)      # Get the neighborhoods of a protein (or of a vector of proteins).
  get_subnetwork(string_ids)     # returns a subgraph from the given input proteins
  
Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        return(graph)
      },
      
      
      
      get_homologs = function(string_ids, target_species_id, bitscore_threshold=NULL){
'
Description:
  Returns the homologs of the given input identifiers that are present in the given target_species_id.

Input parameters:
    "string_ids"          a vector of STRING identifiers.
    "target_species_id"   NCBI taxonomy identifier of the species to query for homologs (the species must be present in the STRING database)
    "bitscore_threshold"  threshold on the bitscore of the blast alignment.

Author(s):
   Andrea Franceschini

'
        
        if(length(string_ids) > 300) {
          cat("ERROR: We support a maximum of 300 STRING identifiers per call. Please reduce the size of the input and try again. \t")
          stop()
        }
        urlStr = paste("http://string-db.org/version_", version, "/newstring_cgi/webservices/homology_hits.pl", sep="")
        identifiers=""
        for(id in string_ids ){
          identifiers = paste(identifiers, id, "%0D", sep="")
        }
        params = list(target_species_id=target_species_id, identifiers=identifiers)
        if(!is.null(bitscore_threshold)) params["bitscore_threshold"] = bitscore_threshold
        tempDfv=postFormSmart(urlStr, .params=params)
        hhits <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        return(hhits)
      },
      
      
      
      
      
      get_homologs_besthits = function(string_ids, symbets = FALSE, target_species_id = NULL, bitscore_threshold=NULL){
'
Description:
  Returns the best blast hits x species of the given input identifiers.

Input parameters:
    "string_ids"            a vector of STRING identifiers.
    "target_species_id"     NCBI taxonomy identifier of the species to query for homologs (the species must be present in the STRING database)
    "bitscore_threshold"    threshold on the bitscore of the blast alignment.
    "symbets"               specify whether you want only symmetrical best hits

Author(s):
   Andrea Franceschini
'
        if(symbets) symbets=1
        else symbets=0
        if(length(string_ids) > 300) {
          cat("ERROR: We support a maximum of 300 STRING identifiers per call. Please reduce the size of the input and try again. \t")
          stop()
        }
        urlStr = paste("http://string-db.org/version_", version, "/newstring_cgi/webservices/homology_best_hits.pl", sep="")
        identifiers=""
        for(id in string_ids ){identifiers = paste(identifiers, id, "%0D", sep="")}
        params = list(identifiers=identifiers, symbets=symbets)
        if(!is.null(bitscore_threshold)) params["bitscore_threshold"] = bitscore_threshold
        if(!is.null(target_species_id)) params["target_species_id"] = target_species_id
        tempDfv=postFormSmart(urlStr, .params=params)
        best_hits <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        best_hits = subset(best_hits, select=c("STRING_id", "species_id", 
                                               "best_hit_STRING_id","best_hit_bitscore", "best_hit_normscore", "best_hit_alignment_length", "nr_high_scoring_hits"))
        return(best_hits)
      },  
      
      
      # find interactions between a group of proteins
      get_interactions = function(string_ids){
'
Description:
  Shows the interactions in between the proteins that are given in input.

Input parameters:
    "string_ids"    a vector of STRING identifiers.

Author(s):
   Andrea Franceschini

'
        if(is.null(graph)) load()
        return(get.data.frame(get_subnetwork(string_ids),  c("edges")))
      },
      
      
      
      get_link = function(string_ids, required_score=NULL, network_flavor="evidence", payload_id = NULL){
'
Description:
  Returns a short link to the network page of our STRING website that shows the protein interactions between the given identifiers.

Input parameters:
  "string_ids"        a vector of STRING identifiers.
  "required_score"    minimum STRING combined score of the interactions 
                        (if left NULL we get the combined score of the object, which is 400 by default)
  "network_flavor"    specify the flavor of the network ("evidence", "confidence" or "actions".  default "evidence")

Author(s):
   Andrea Franceschini
'
        if(length(string_ids) > 400) {
          cat("ERROR: We do not support lists with more than 400 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
          stop()
        }
        if(is.null(required_score) ) required_score = score_threshold
        string_ids = unique(string_ids)
        urlStr = paste("http://string-db.org/version_",version,"/newstring_cgi/link_to.pl", sep="" )
        #urlStr = paste("http://130.60.240.82:8818/newstring_cgi/link_to.pl", sep="" )
        identifiers=""
        for(id in string_ids ){   identifiers = paste(identifiers, id, "%0D", sep="")}
        params=list(required_score=required_score, limit=0, network_flavor=network_flavor, identifiers=identifiers)
        if(!is.null(payload_id)) params["internal_payload_id"] = payload_id
        tempDfv=postFormSmart(urlStr, .params=params)
        df=read.table(text=tempDfv, stringsAsFactors=FALSE, fill = TRUE)
        return(df$V1)

#           if(is.null(required_score) ) required_score = score_threshold
#           string_ids = unique(string_ids)
#           urlStr = paste("http://string-db.org/version_",version,"/newstring_cgi/link_to.pl?required_score=", required_score, "&limit=0&network_flavor=", network_flavor, sep="" )
#           if(!is.null(payload_id)) urlStr = paste(urlStr, "&internal_payload_id=", payload_id, sep="")
#           urlStr = paste(urlStr, "&identifiers=", sep="")
#           for(id in string_ids ){
#             urlStr = paste(urlStr, id, "%0D", sep="")
#           }
#           df=read.table(url(urlStr), stringsAsFactors=FALSE, fill = TRUE)
#           return(df$V1)


      },
      
      
      # find interactors of a specific protein
      get_neighbors = function(string_ids){
'
Description:
Get the neighborhoods of a protein (or of a vector of proteins) that is given in input.

Input parameters:
  "string_ids" =  a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        return(V(graph)[neighbors(graph, string_ids)]$name)
      },
      
      
      
      
      get_png = function(string_ids, required_score=NULL, network_flavor="evidence", file=NULL, payload_id=NULL){
'
Description:
  Returns a png image of a STRING protein network with the given identifiers.

Input parameters:
  "string_ids"        a vector of STRING identifiers.
  "required_score"    minimum STRING combined score of the interactions 
                        (if left NULL we get the combined score of the object, which is 400 by default)
  "network_flavor"    specify the flavor of the network ("evidence", "confidence" or "actions".  default "evidence")
  "file"              file where to save the image (must have .png extension)

Author(s):
   Andrea Franceschini
'
        if(length(string_ids) > 400) {
          cat("ERROR: We do not support lists with more than 400 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
          stop()
        }
        if(is.null(required_score) ) required_score = score_threshold
        string_ids = unique(string_ids)
        string_ids = string_ids[!is.na(string_ids)]
        urlStr = paste("http://string-db.org/version_", version, "/api/image/network", sep="" )
        identifiers=""
        for(id in string_ids ){ identifiers = paste(identifiers, id, "%0D", sep="")}
        params = list(output="image", required_score=required_score, limit=0, network_flavor=network_flavor, identifiers=identifiers)
        if(!is.null(payload_id)) params["internal_payload_id"]= payload_id
        img <- readPNG(postFormSmart(  urlStr, .params=params) )
        if(!is.null(file))  writePNG(img,  file)
        
        return(img)
      },
      
      
      
      get_proteins = function(){
'
Description:
  Returns the STRING proteins data frame.
  (it downloads and caches the information the first time that is called).

Author(s):
   Andrea Franceschini
'
        
        if(nrow(proteins)==0){
          temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/proteins/",
                                          species, "__proteins.tsv.gz", sep=""), oD=input_directory)
          proteinsDf <- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE, quote="")
          proteinsDf2 = subset(proteinsDf, select=c("protein_external_id",  "preferred_name", "protein_size", "annotation"))
          proteins <<- proteinsDf2
        }
        return(proteins)
      },
      
      
      
      
      get_ppi_enrichment = function(string_ids){
'
Description:
  Returns a pvalue representing the enrichment in interactions of the list of proteins 
          (i.e. the probability to obtain such a number of interactions by chance).

Input parameters:
  "string_ids"    a vector of STRING identifiers, sorted.

Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        return(ppi_enrichment(string_ids, graph))
      },
      
      
      get_ppi_enrichment_full = function(string_ids, sliceWindow = 20, edgeWindow  = 140, windowExtendedReferenceThreshold = 260, growingWindowLimit=NULL){
'
Description:
  Returns a vector showing the enrichment in protein interactions in various positions of the list of genes in input.
  In practice, a list of 3 vectors is returned:
  1) enrichment  (i.e.  enrichment computed in the window from 1 to x)
  2) enrichmentWindow (i.e. enrichment computed in a sliding window of size determined by the "edgeWindow" parameters 
  and the sliding steps determined by the "sliceWindow" parameter)
  3) enrichmentWindowExtended  (i.e. like the enrichmentWindow, 
  but it also includes an initial window of size "windowExtendedReferenceThreshold" with respect to which to compute the enrichment )

Input parameters:
  "string_ids"                          vector of STRING identifiers, sorted.
  "file"                                file where to save the graph as an image
  "sliceWindow"                         defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
  "windowExtendedReferenceThreshold"    defines the size of a window at the beginning of the list. 
                                              The enrichment will be computed always including the proteins in this window
  "title"                               title of the graph.
  "growingWindowLimit"                  threshold where to stop the computation of the enrichment

Author(s):
   Andrea Franceschini

'  
        if(is.null(graph)) load()
        return(ppi_enrichment_full(string_ids, graph, sliceWindow = sliceWindow, edgeWindow  = edgeWindow, 
                                   windowExtendedReferenceThreshold = windowExtendedReferenceThreshold, growingWindowLimit=growingWindowLimit))
      },
      
      
#       get_species = function(species_name=NULL){
# '
# Input parameters:
#   "species_name"      name of the species that you are searching (e.g. "Homo")
# Returns a list of species that are present in STRING.
# '
#         
#         if(nrow(speciesList)==0){
#           temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/species.tsv", sep=""), oD=input_directory)
#           speciesList <<- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)  
#         }
#         speciesDf = speciesList
#         if(!is.null(species_name)){
#           speciesDf = subset(speciesDf, grepl(species_name, official_name))
#         }
#         return(speciesDf)
#       },
      
      
      get_pubmed = function(string_ids){
'
Description:
  Returns vector with the PUBMED IDs of the publications that contain the names of the proteins in the input vector.

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'        
        if(length(string_ids) > 300) {
          cat("ERROR: We support a maximum of 300 STRING identifiers per call. Please reduce the size of the input and try again. \t")
          stop()
        }
        #urlStr = paste("http://string-db.org/version_", version, "/api/tsv/abstracts?limit=1000000&identifiers=", sep="")
        urlStr = paste("http://string-db.org/version_", version, "/newstring_cgi/webservice_handler.pl", sep="")
        identifiers=""
        for(id in string_ids ){identifiers = paste(identifiers, id, "%0D", sep="")}
        params = list(limit=1000000, identifiers=identifiers, output="tsv", request="abstracts")
        tempDfv=postFormSmart(urlStr, .params=params)
        pubmedIdsDf <- read.table(text=tempDfv, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        return(pubmedIdsDf$abstractId)
      },
      
      
      get_pubmed_interaction = function(STRING_id_a, STRING_id_b){
'
Description:
  Returns vector with the PUBMED IDs of the publications that contain the names of both the input proteins.

Input parameters:
  "STRING_id_a"      STRING identifier.
  "STRING_id_b"      STRING identifier.

Author(s):
   Andrea Franceschini
'             
        return(intersect(get_pubmed(STRING_id_a), get_pubmed(STRING_id_b) ) )
      },
      
      
      get_subnetwork = function(string_ids){
'
Description:
  Returns the subgraph generated by the given input proteins.

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        return(induced.subgraph(graph, which(V(graph)$name %in% string_ids)))
      },
      
      
      get_summary = function(string_ids){
'
Description:
  Returns a summary of the STRING sub-network containing the identifiers provided in input.

Input parameters:
  "string_ids"      a vector of STRING identifiers.

Author(s):
   Andrea Franceschini
'
        if(is.null(graph)) load()
        hitsWithEdges = string_ids[string_ids %in% V(graph)$name]
        enrichment = ppi_enrichment(string_ids, graph)
        summaryText = paste("proteins: ", length(string_ids), '\n', 
                            "interactions: ", length(E(induced.subgraph(graph, hitsWithEdges ))),'\n',
                            "expected interactions: ", enrichment$lambda,' (',
                            "p-value: ",enrichment$enrichment, ')\n', sep=""
        )
        return(summaryText)
      },
      
      
      
      get_term_proteins = function(term_ids, string_ids=NULL, enableIEA = TRUE){
'
Description:
  Returns the proteins annotated to belong to a given term.

Input parameters:
  "term_ids"        vector of terms 
  "string_ids"      vector of STRING identifiers. If the variable is set, the method returns only the proteins that are present in this vector.
  "enableIEA"       whether to consider also Electronic Inferred Annotations

Author(s):
   Andrea Franceschini
'
        
        annotations2 = get_annotations()
        if(!enableIEA) { annotations2 = subset(annotations2, type!="IEA") }
        annotations3 = subset(annotations2, term_id %in% term_ids, select=c("STRING_id", "term_id"))
        if(!is.null(string_ids)){
          annotations3 = subset(annotations3, STRING_id %in% string_ids)
        }
        annotations4 = add_proteins_description(annotations3)
        return(annotations4)
      },
      
      
      load = function(){
'
Description:
  Downloads and returns the STRING network (the network is set also in the graph variable of the STRING_db object).
  
It makes use of the variables:
    "backgroundV"         vector containing STRING identifiers to be used as background 
                            (i.e. the STRING network loaded will contain only the proteins that are present also in this vector)
    "score_threshold"     STRING combined score threshold (the network loaded contains only interactions having a combined score greater than this threshold)

Author(s):
   Andrea Franceschini
'
        
        
        temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/protein_links/",
                                        species, "__protein_links.tsv.gz", sep=""), oD=input_directory)
        PPI <- read.table(temp, sep = " ", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)
        
        PPIselected = PPI
        if(length(score_threshold)!=0) PPIselected <- PPI[PPI$combined_score >= score_threshold,]
        if(!is.null(backgroundV)){
          PPIselected <- PPIselected[PPIselected$protein1 %in% backgroundV,]
          PPIselected <- PPIselected[PPIselected$protein2 %in% backgroundV,]
        }
        
        myg = graph.data.frame(PPIselected,FALSE)
        graph <<- myg
        return(myg)
      },
      
      
      
      mp = function(protein_aliases){
'
Description:
Maps the gene identifiers of the input vector to STRING identifiers (using a take first approach).
It returns a vector with the STRING identifiers of the mapped proteins

Input parameters:
    "protein_aliases"       vector with the aliases of the proteins that we want to map to STRING

Author(s):
Andrea Franceschini
'
        
        temp_df = data.frame(proteins=protein_aliases)
        temp_df_mapped = map(temp_df, "proteins", removeUnmappedRows=TRUE, quiet=TRUE)
        return(temp_df_mapped$STRING_id)
        
      },
      
      
      map = function(my_data_frame,
                              my_data_frame_id_col_names,
                              takeFirst=TRUE, removeUnmappedRows=FALSE, quiet=FALSE
      ){
'
Description:
  Maps the gene identifiers of the input dataframe to STRING identifiers.
  It returns the input dataframe with the "STRING_id" additional column.

Input parameters:
  "my_data_frame"                 data frame provided as input. 
  "my_data_frame_id_col_names"    vector contatining the names of the columns of "my_data_frame" that have to be used for the mapping.
  "takeFirst"                     boolean indicating what to do in case of multiple STRING proteins that map to the same name. 
                                      If TRUE, only the first of those is taken. Otherwise all of them are used. (default TRUE)
  "removeUnmappedRows"            remove the rows that cannot be mapped to STRING 
                                      (by default those lines are left and their STRING_id is set to NA)
  "quiet"                         Setting this variable to TRUE we can avoid printing the warning relative to the unmapped values.

Author(s):
   Andrea Franceschini
'

        aliasDf2=NULL
        if(!takeFirst || length(aliases_tf)==0){        
          if(takeFirst)
            temp = downloadAbsentFileSTRING(paste("http://string.uzh.ch/permanent/string/", version, "/protein_aliases/", species, "__protein_aliases_tf.tsv.gz", sep=""), oD=input_directory)
          else
            temp = downloadAbsentFileSTRING(paste("http://string.uzh.ch/permanent/string/", version, "/protein_aliases/", species, "__protein_aliases.tsv.gz", sep=""), oD=input_directory)
          
          aliasDf <- read.table(temp, sep = "\t", header=TRUE, quote="", stringsAsFactors=FALSE, fill = TRUE)
          aliasDf = renameColDf(aliasDf, "protein_id", "STRING_id")
          aliasDf2 = subset(aliasDf, select=c("STRING_id", "alias"))
          if(takeFirst) aliases_tf <<- aliasDf2
        }else{
          aliasDf2 = aliases_tf
        }


        tempDf = multi_map_df(my_data_frame, aliasDf2, my_data_frame_id_col_names, "alias", "STRING_id")
        naDf = subset(tempDf, is.na(STRING_id))
        if(nrow(naDf) > 0 & !quiet) cat(paste("Warning:  we couldn't map to STRING ", as.integer((nrow(naDf)/nrow(tempDf))*100), "% of your identifiers" , sep=""))
        if(removeUnmappedRows) tempDf = subset(tempDf, !is.na(STRING_id))
        return(tempDf)
        
      },
      
      
      plot_ppi_enrichment = function(string_ids, file=NULL, sliceWindow = 20, edgeWindow = 140, 
                                     windowExtendedReferenceThreshold = 260, minVal=0.0000000001, title="", quiet=FALSE){
'
Description:
  Plots a graph showing the enrichment in protein interactions in various positions of the list of genes in input.

Input parameters:
  "string_ids"                          vector of STRING identifiers, sorted.
  "file"                                file where to save the graph as an image
  "sliceWindow"                         defines the interval in proteins after which to compute the enrichment, scanning the list (i.e. the resolution)
  "windowExtendedReferenceThreshold"    defines the size of a window at the beginning of the list. 
                                          The enrichment will be computed always including the proteins in this window.
  "minVal"                              minimum value of the pvalue (lower values with respect to this one will assume this minimum value.)
  "title"                               title of the graph.

Author(s):
   Andrea Franceschini

'
          if(is.null(graph)) load()
          plot_ppi_enrichment_graph(string_ids, graph, file=file, sliceWindow = sliceWindow, edgeWindow = edgeWindow, 
                                    windowExtendedReferenceThreshold = windowExtendedReferenceThreshold, minVal=minVal, title=title, quiet=quiet)
      },
      
      
      plot_network = function(string_ids, payload_id=NULL, required_score=NULL, add_link=TRUE, add_summary=TRUE){
'
Description:
  Plots an image of the STRING network with the given proteins.

Input parameters:
  "string_ids"        a vector of STRING identifiers
  "payload_id"        an identifier of payload data on the STRING server (see method post_payload for additional informations)
  "score_threshold"   a threshold on the score that overrides the default score_threshold, that we use only for the picture
  "add_link"          parameter to specify whether you want to generate and add a short link to the relative page in STRING. 
                      As default this option is active but we suggest to deactivate it in case one is generating many images (e.g. in a loop). 
                      Deactivating this option avoids to generate and store a lot of short-urls on our server.
  "add_summary"       parameter to specify whether you want to add a summary text to the picture. This summary includes a p-value and the number of proteins/interactions.

Author(s):
   Andrea Franceschini
'
        
        if(is.null(required_score) ) required_score = score_threshold
        img = get_png(string_ids, payload_id=payload_id, required_score=required_score)
        if(!is.null(img)){
          plot(1:(dim(img)[2]), type='n', xaxt='n', yaxt='n', xlab="", ylab="", ylim=c(1,dim(img)[1]), xlim=c(1,(dim(img)[2])), asp = 1 )
          if(add_summary) mtext(get_summary(string_ids), cex = 0.4)
          if(add_link) mtext(get_link(string_ids, payload_id=payload_id, required_score=required_score), cex = 0.4, side=1)
          rasterImage(img, 1, 1, (dim(img)[2]), dim(img)[1])
        } 
      },
      
      
      post_payload = function(stringIds, colors=NULL, comments=NULL, links=NULL, iframe_urls=NULL, logo_imgF=NULL, legend_imgF=NULL ){
'
Description:
  Posts the input to STRING and returns an identifier that you can use to access the payload when you enter in our website.

Input parameters:
  "string_ids"        vector of STRING identifiers.
  "colors"            vector containing the colors to use for a every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "comments"          vector containing the comments to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "links"             vector containing the links to use for every STRING identifier ( the order of the elements must match those in the string_ids vector)
  "iframe_urls"       vector containing the urls of the iframes to use for every STRING identifier 
                              ( the order of the elements must match those in the string_ids vector)
  "logo_imgF"         path to a file containing the logo image to be display in the STRING website
  "legend_imgF"       path to a file containing a legend image to be display in the STRING website

Author(s):
   Andrea Franceschini

'
        
        postFormParams = list(identifiers=paste(stringIds, collapse=" ") )
        if(!is.null(colors)) postFormParams = c(postFormParams, colors=paste(colors, collapse=" "))
        if(!is.null(comments)) postFormParams = c(postFormParams, comments=paste(comments, collapse=" "))
        if(!is.null(links)) postFormParams = c(postFormParams, links=paste(links, collapse=" "))
        if(!is.null(iframe_urls)) postFormParams = c(postFormParams, iframe_urls=paste(iframe_urls, collapse=" "))
        if(!is.null(logo_imgF)) postFormParams = c(postFormParams, logo_img=fileUpload(logo_imgF))
        if(!is.null(legend_imgF)) postFormParams = c(postFormParams, legend_img=fileUpload(legend_imgF))
        postRs = postFormSmart(paste("http://string-db.org/version_", version, "/newstring_cgi/webservices/post_payload.pl", sep=""),  .params = postFormParams)
        return(postRs)
      },
      
      
      set_background = function(background_vector){
'
Description:
  With this method you can specify a vector of proteins to be used as background. 
  The network is reloaded and only the proteins that are present in the background vector are inserted in the graph.  
  Besides, the background is taken in consideration for all the enrichment statistics.

Input parameters:
  "background_vector"     vector of STRING protein identifiers

Author(s):
   Andrea Franceschini
'
        
        backgroundV <<- background_vector
        load()
      },
      
  
      show = function(){
        cat(paste(
          "***********  STRING - http://string-db.org   ***********", "\n",
          "(Search Tool for the Retrieval of Interacting Genes/Proteins)  ", "\n",
          "version: ", version, "\n",
          "species: ", species, "    ", subset(get_STRING_species(), species_id==species)$official_name, "\n", 
          "............please wait............\n", sep=""
        ))
        
        cat(paste(
            "proteins: ", nrow(get_proteins()), "\n",
            "interactions: ", length(E(get_graph())), sep=""
                )
            )
        
      }
    
      
    )                     
)



get_STRING_species = function(version="9_05", species_name=NULL){
  
#   Input parameters:
#   "species_name"      name of the species that you are searching (e.g. "Homo")
#   
#   Description:
#   Returns a list of species that are present in STRING.
# 
#   Author(s):
#   Andrea Franceschini
  
  
  
  temp = downloadAbsentFile(paste("http://string.uzh.ch/permanent/string/", version, "/species.tsv", sep=""), oD=tempdir())
  speciesList <- read.table(temp, sep = "\t", header=TRUE, stringsAsFactors=FALSE, fill = TRUE)  
  speciesDf = speciesList
  if(!is.null(species_name)){
    speciesDf = subset(speciesDf, grepl(species_name, official_name))
  }
  return(speciesDf)
}






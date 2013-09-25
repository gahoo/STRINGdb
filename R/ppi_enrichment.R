#
#  copyright:   Andrea Franceschini
#          (Swiss Institute of Bioinformatics)
#           andrea.franceschini@isb-sib.ch
#



library(igraph)


ppi_enrichment <- function( hitList, ppi_network){
  hitListIndices <- rep(NA, length(hitList))
  for( i in 1:length(hitList) ){ if(hitList[i] %in% V(ppi_network)$name) {hitListIndices[i] <- which(V(ppi_network)$name==hitList[i]);} else {hitListIndices[i] <- (-1)} }
  hitListSlice = hitListIndices
  hitListSliceEdgeNumber = ecount(induced.subgraph(ppi_network, hitListSlice[hitListSlice>=0]))
  numNodes = length(hitListSlice[hitListSlice>=0])
  totalEdgeNum = ecount(ppi_network)
  lambda = ppie.compLambda(degree(ppi_network, hitListSlice[hitListSlice>=0]), totalEdgeNum)
  enrichment = ( ppois(lambda=lambda, ((numNodes*(numNodes-1))/2), lower.tail=TRUE ) - 
    ppois(lambda=lambda, (hitListSliceEdgeNumber-1), lower.tail=TRUE ) ) / 
    ppois(lambda=lambda, ((numNodes*(numNodes-1))/2), lower.tail=TRUE )  
  
  tempReturn = list(
    enrichment = enrichment,
    lambda = as.integer(lambda)
  )
  
  return(tempReturn)
}

#
# Computes the enrichment in protein-protein interactions of a list of proteins.
# 
# Input parameters:
#   "hitList" = List of proteins for which to compute the enrichment, in various positions of the list
#   "ppi_network" = protein-protein interactions graph
#   "sliceWindow" = size of the window that we use to scan the input list
#               (i.e. a pvalue is computed each "sliceWindow" proteins. )
#   "edgeWindow" = size of the window that we use to compute the enrichment (i.e. the window pvalue is computed using the proteins inside this "edgeWindow")
#   "windowExtendedReferenceThreshold" = when we compute the "windowExtended" pvalue we are computing the pvalue that considers the following interactions:
#                                             1) the interactions inside the edgeWindow (as we do with the edgeWindow pvalue)
#                                             2) the interactions that connects the proteins in the edgeWindow with the proteins in another window at the beginning of the list (i.e. the windowExtendedReference).
#                                         windowExtendedReferenceThrehold defines the size of this windowExtendedReference window.
#                                         In this way we can compute, in a reliable way, the enrichment of a sorted list of proteins, in various positions of the list.
#   "growingWindowLimit" = stop to compute the enrichment (from position 1 to position n) after growingWindowLimit proteins in the sorted list.
#                                          (this limit speeds up the computation of the 2 other types of enrichment)
ppi_enrichment_full <- function( hitList, ppi_network, sliceWindow = 20, edgeWindow  = 140, windowExtendedReferenceThreshold = 260, growingWindowLimit=NULL, quiet=FALSE){
  slices=as.integer(length(hitList)/sliceWindow)
  if(!quiet) pb <- txtProgressBar(min = 0, max = slices, style = 3)
  hitListIndices <- rep(NA, length(hitList))
  for( i in 1:length(hitList) ){ if(hitList[i] %in% V(ppi_network)$name) {hitListIndices[i] <- which(V(ppi_network)$name==hitList[i]);} else {hitListIndices[i] <- (-1)} }
  
  enrichment <- rep(NA, slices)
  enrichmentWindow <- rep(NA, slices)
  enrichmentWindowExtended <- rep(NA, slices)
  
  totalEdgeNum = ecount(ppi_network)
  
  for(i in 1:slices){
    if(!quiet) setTxtProgressBar(pb, i)
    limit=sliceWindow*i
    
    if(i<=edgeWindow/sliceWindow || ((!is.null(growingWindowLimit)) && i<growingWindowLimit) ){
      hitListSlice = hitListIndices[1:limit]
      hitListSliceEdgeNumber = ecount(induced.subgraph(ppi_network, hitListSlice[hitListSlice>=0]))
      numNodes = length(hitListSlice[hitListSlice>=0])
      lambda = ppie.compLambda(degree(ppi_network, hitListSlice[hitListSlice>=0]), totalEdgeNum)
      
      enrichment[i] = ( ppois(lambda=lambda, ((numNodes*(numNodes-1))/2), lower.tail=TRUE ) - 
        ppois(lambda=lambda, (hitListSliceEdgeNumber-1), lower.tail=TRUE ) ) / 
        ppois(lambda=lambda, ((numNodes*(numNodes-1))/2), lower.tail=TRUE )
    }
    
    if(i>edgeWindow/sliceWindow) {
      limitBack = sliceWindow*i-edgeWindow
      hitListSliceBack = hitListIndices[1: min(limitBack, windowExtendedReferenceThreshold)]
      hitListSliceWindow = hitListIndices[(limitBack+1):limit]
      hitListSliceBack = hitListSliceBack[hitListSliceBack>=0]
      hitListSliceWindow = hitListSliceWindow[hitListSliceWindow>=0]
      numNodesW = length(hitListSliceWindow)
      numNodesBack = length(hitListSliceBack)
      numEdgesW = ecount(induced.subgraph(ppi_network, hitListSliceWindow))
      
      hitListSliceWindowDegrees = degree(ppi_network, hitListSliceWindow)
      lambdaInsideW = ppie.compLambda(hitListSliceWindowDegrees, totalEdgeNum)
      lambdaBeforeW = ppie.compLambdaL1L2(degree(ppi_network, hitListSliceBack), hitListSliceWindowDegrees, totalEdgeNum);
      lambdaW = lambdaInsideW + lambdaBeforeW;
      
      hitListSliceEdgesNumW = numEdgesW + ppie.getNumEdgesBetween(ppi_network, hitListSliceWindow, hitListSliceBack)
      maxNumEdgesExtendedW = ((numNodesW*(numNodesW-1))/2) + numNodesW*numNodesBack   
      
      if(limit < windowExtendedReferenceThreshold){
        hitListSliceForward = hitListIndices[(limit+1):windowExtendedReferenceThreshold]
        hitListSliceForward = hitListSliceForward[hitListSliceForward>=0]
        hitListSliceForwardDegrees = degree(ppi_network, hitListSliceForward)
        lambdaAfterW = ppie.compLambdaL1L2(hitListSliceForwardDegrees, hitListSliceWindowDegrees, totalEdgeNum)
        lambdaW = lambdaW + lambdaAfterW
        
        hitListSliceEdgesNumW = hitListSliceEdgesNumW + ppie.getNumEdgesBetween(ppi_network, hitListSliceWindow, hitListSliceForward)
        maxNumEdgesExtendedW = maxNumEdgesExtendedW + ( length(hitListSliceForwardDegrees) * numNodesW )
      }
      
      enrichmentWindowExtended[i] = ( ppois(lambda=lambdaW, (maxNumEdgesExtendedW), lower.tail=TRUE ) - 
        ppois(lambda=lambdaW, (hitListSliceEdgesNumW-1), lower.tail=TRUE ) ) / 
        ppois(lambda=lambdaW, (maxNumEdgesExtendedW), lower.tail=TRUE )
      
      maxNumEdgesW = ( numNodesW * (numNodesW-1) ) / 2;
      enrichmentWindow[i] = ( ppois(lambda=lambdaInsideW, (maxNumEdgesW), lower.tail=TRUE ) - 
        ppois(lambda=lambdaInsideW, (numEdgesW-1), lower.tail=TRUE ) ) / 
        ppois(lambda=lambdaInsideW, (maxNumEdgesW), lower.tail=TRUE )
      
    }else{
      enrichmentWindow[i]=NA
      enrichmentWindowExtended[i] = enrichment[i]
    }
    
  }
  
  tempReturn = list(
    enrichment = enrichment,
    enrichmentWindow = enrichmentWindow,
    enrichmentWindowExtended = enrichmentWindowExtended
  )
  
  if(!quiet) close(pb)
  return(tempReturn)
}

#
# Plots a graph with the enrichment in protein interactions
# See ppi_enrichment_full comments for additional details.
#
plot_ppi_enrichment_graph <- function(proteins, ppi_network, file=NULL, sliceWindow = 20, edgeWindow = 140, 
                                          windowExtendedReferenceThreshold = 260, minVal=0.0000000001, title="", quiet=FALSE){
  
  enrichment = ppi_enrichment_full(proteins, ppi_network, sliceWindow = sliceWindow, edgeWindow  = edgeWindow, windowExtendedReferenceThreshold = windowExtendedReferenceThreshold, quiet=quiet)
  enrichment2 = pmin(pmax(rep(minVal, length(enrichment$enrichmentWindowExtended)), enrichment$enrichmentWindowExtended*15), 1)
  
  if(!is.null(file))
    if(regexpr("pdf", file)>0) pdf(file)
  else if(regexpr("jpg", file)>0) jpeg(file)
  
  if(title=="") title = "PPI Enrichment"
  
  plot(y=log(enrichment2[0:length(enrichment2)], base=10), x=seq(0, sliceWindow*(length(enrichment2)-1), 20), type="b", col="blue", xlab="index", 
       ylab="log(p-value)", ylim=c(-10,0), main=title )
  abline(h=log(0.05, base=10))
  #par(xpd=T)
  #legend(-200, 100, c("Enrichment"), col=c("blue"), pch=c(21))
  #par(xpd=F)
  if(!is.null(file)) dev.off()
  
}


ppie.getNumEdgesBetween <- function(graph, nodesFrom, nodesTo){
  edgesEnd = unlist(neighborhood(graph, 1, nodesFrom))
  return(length( edgesEnd[edgesEnd %in% nodesTo] ) )
}


# ppie.getNeighborhoods <- function(graph, inputNodes){
#   return(unique(unlist(neighborhood(graph, 1, inputNodes))))
# }


ppie.compLambdaL1L2 <- function(degreesI, degreesJ, edgeNum){
  lambda = 0
  for(i in 1:length(degreesI)){
    for(j in 1:length(degreesJ)){
      lambda = lambda + ppie.compPij(degreesI[i], degreesJ[j], edgeNum)
    }
  }
  return(lambda)
}


ppie.compLambda <- function(degrees, edgeNum){
  lambda = 0
  for(i in 1:length(degrees)){
    for(j in 1:i){
      lambda = lambda + ppie.compPij(degrees[i], degrees[j], edgeNum)
    }
  }
  return(lambda)
}


ppie.compPij <- function(degI, degJ, edgeNum){
  compPijVal = 1-exp(-((degI*degJ)/(2*edgeNum)))
  return(compPijVal)
}


#
#  copyright:   Andrea Franceschini
#          (Swiss Institute of Bioinformatics)
#           andrea.franceschini@isb-sib.ch
#


postFormSmart <- function(uri, ..., .params = list(), .opts = curlOptions(url = uri),
                          curl = getCurlHandle(), style = 'HTTPPOST',
                          .encoding = integer(), binary = NA, .checkParams = TRUE,
                          .contentEncodeFun = curlEscape){

    res = postForm(uri, ..., .params = .params, .opts = .opts,
             curl = curl, style = style,
             .encoding = .encoding, binary = binary, .checkParams = .checkParams,
             .contentEncodeFun = .contentEncodeFun)
  
    
    suppressWarnings( if(grepl("The document has moved", res)){
    
      begin <- regexpr("href",res)+6
      mys2=substr(res, begin, 10000000)
      end <- regexpr('"',mys2)-1
      uriNew = substr(mys2, 1, end)
      
      res=postForm(uriNew, ..., .params = .params, .opts = .opts,
               curl = curl, style = style,
               .encoding = .encoding, binary = binary, .checkParams = .checkParams,
               .contentEncodeFun = .contentEncodeFun)
    } )
    
    return(res)
  
}


coeffOfvar <- function(x){
  return(sd(x)/mean(x))
} 


# delete column in data frame
delColDf <- function(df, colName){
  if(colName %in% names(df)) return( df[ , -which(names(df) %in% c(colName))] )
  else return(df)
}


# download a file when it is not already present and also check the dimension of the file with the STRING server
downloadAbsentFileSTRING <- function(urlStr, oD = tempdir()){
  expectedSize=read.table(url(paste("http://string.uzh.ch/permanent_scripts/get_file_size.pl?file=", urlStr, sep="")), stringsAsFactors=FALSE, fill = TRUE)$V1
  
  fileName = tail(strsplit(urlStr, "/")[[1]], 1)
  temp <- paste(oD,"/", fileName, sep="")
  if(! file.exists(temp) || file.info(temp)$size==0 || file.info(temp)$size!=expectedSize) download.file(urlStr,temp)
  if(file.info(temp)$size==0) {
    unlink(temp)
    temp=NULL
    cat(paste("ERROR: failed to download ", fileName,".\nPlease check your internet connection and/or try again. " , 
              "\nThen, if you still display this error message please contact us.",sep=""))
  }  
  return(temp)
}


# download a file when it is not already present
downloadAbsentFile <- function(urlStr, oD = tempdir()){
  
  fileName = tail(strsplit(urlStr, "/")[[1]], 1)
  temp <- paste(oD,"/", fileName, sep="")
  if(! file.exists(temp) || file.info(temp)$size==0) download.file(urlStr,temp)
  if(file.info(temp)$size==0) {
    unlink(temp)
    temp=NULL
    cat(paste("ERROR: failed to download ", fileName,".\nPlease check your internet connection and/or try again. " , 
              "\nThen, if you still display this error message please contact us.",sep=""))
  }  
  return(temp)
}



# drop levels from vector
# dropLevelsFV <- function(vect){
#   newVect = rep(NA, length(vect))
#   for(i in 1:length(vect)){newVect[i] = vect[i]}
#   return(newVect)
# } 


# sort the vector and extract the values in between startPerc and endPerc 
# extract_values <- function(vect, startPerc, endPerc){
# 	vect = vect[!is.na(vect)]
# 	indexA=floor((length(vect)* startPerc) + 1)
# 	indexB=floor(length(vect)* endPerc)
# 	vectToReturn = sort(vect)[indexA:indexB]
# 	return(vectToReturn)
# }


# intersectAll <-function(...){
#   vectors <- list(...)
#   intersection = vectors[[1]]
#   for(i in 1:length(vectors)){
#     intersection = intersect(intersection, vectors[[i]])
#   }
#   return(intersection)
# }


## move column in the defined position
# moveCol <- function(df, colName, new_position=1){
# 	col_idx <- grep(colName, names(df))
# 	posArr = 1:ncol(df)
# 	posArr = append(posArr[-col_idx], col_idx, after=new_position-1)
# 	df <- df[, posArr]
# 	return(df)
# }


## move column in the position after colNameBefore
# moveColAfter <- function(df, colName, colNameBefore){
# 	new_position <- (grep(colNameBefore, names(df)))+1
# 	return(moveCol(df,colName,new_position))
# }



# Merge preserving the original order if sort==FALSE (this doesn't happen in the original merge implementation)
merge.with.order <- function(x,y, ..., sort = T)
{
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }
  
  if(sort==F){ return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
  } else {return(merge(x=x,y=y,..., sort = sort))}
  
}



# mapping function (add the possibility to map using more than one column)
multi_map_df <- function(dfToMap, dfMap, strColsFrom, strColFromDfMap, strColToDfMap, caseSensitive=FALSE){
  
  tempMatr = matrix(NA, length(strColsFrom), nrow(dfToMap))
  for(i in 1:length(strColsFrom)){
    if(!caseSensitive){ 
      tempMatr[i,] = as.vector(dfToMap[,strColsFrom[i]])
      dfToMap[,strColsFrom[i]] = toupper(iconv(dfToMap[,strColsFrom[i]], "WINDOWS-1252", "UTF-8")) 
    }
  }
  if(!caseSensitive){ dfMap[,strColFromDfMap] = toupper(iconv(dfMap[,strColFromDfMap], "WINDOWS-1252", "UTF-8")) }
  
  dfMap2 = unique(subset(dfMap, select=c(strColFromDfMap, strColToDfMap)))
  df2 = merge.with.order(dfToMap, dfMap2, by.x=strColsFrom[1], by.y=strColFromDfMap, all.x=TRUE, sort=FALSE)
  if(length(strColsFrom) > 1){
    for(i in 2:length(strColsFrom)){
      dfna = delColDf(subset(df2, is.na(as.vector(df2[, strColToDfMap]))), strColToDfMap)
      dfgood = subset(df2, !is.na(as.vector(df2[, strColToDfMap])))
      df3 = merge.with.order(dfna, dfMap2, by.x=strColsFrom[i], by.y=strColFromDfMap, all.x=TRUE, sort=FALSE)
      df2 = rbind(dfgood, df3)
    }  
  }
  
  for(i in 1:length(strColsFrom)){
    if(!caseSensitive && length(tempMatr[i,])==length(df2[,strColsFrom[i]])) df2[,strColsFrom[i]] = tempMatr[i,]
  }
  
  return(df2)
}



# sum n random numbers
# sumRandom <- function(n){
#   rand=0
#   for(i in 1:n){
#     rand=rand+sample(1:100,1,replace=TRUE)/100
#   }
#   return(rand)
# }

  
# generate a vector contaning as elements sums of n random numbers 
# sumRandomVect <- function(n, vLength){
#   vRand = NULL
#   for(i in 1:vLength){
#     vRand = c(vRand, sumRandom(n))
#   }
#   return(vRand)
# }


## Rename column data frame
renameColDf <- function(df, colOldName, colNewName){
	if(! (colOldName %in% names(df)) ) print(paste("ERROR: We cannot find ", colOldName, " in the data frame.", sep=""))
  names(df)[ which( names( df ) == colOldName ) ] <- colNewName
	return(df)
}


# replace elements of inputVect with the corresponding element in replacementVect when it is not null
# replace_non_null_elements <- function(inputVect, replacementVect){
#   for(i in 1:length(inputVect)){
#     if(!is.na(replacementVect[i]) && !is.null(replacementVect[i])) inputVect[i] = replacementVect[i]
#   }
#   return(inputVect)
# }


## round values of the specified columns
# roundDf <- function(df, colNames, digits=3){
# 	for( i in 1:length(colNames) ){
# 		df[colNames[i]] = round(df[colNames[i]], digits)
# 	}
# 	return(df)
# } 



# Union all the vectors
# unionAll <-function(...){
#   vectors <- list(...)
#   unionVectors = vectors[[1]]
#   for(i in 1:length(vectors)){
#     unionVectors = union(unionVectors, vectors[[i]])
#   }
#   return(unionVectors)
# }



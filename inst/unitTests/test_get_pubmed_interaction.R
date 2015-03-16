test_get_pubmed_interaction <- function(){
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=400 ) 
  pi = string_db$get_pubmed_interaction(string_db$mp("tp53"), string_db$mp("atm"))
  checkTrue(!is.null(pi) && length(pi) > 0)
}
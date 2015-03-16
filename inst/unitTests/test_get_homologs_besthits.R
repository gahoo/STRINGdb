test_get_homologs_besthits <- function() {
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=400 ) 
  data(diff_exp_example1)
  example1_mapped = string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )  
  homologs = string_db$get_homologs_besthits( subset(example1_mapped, gene=="TP53")$STRING_id )
  checkTrue(nrow(homologs) > 0)
}
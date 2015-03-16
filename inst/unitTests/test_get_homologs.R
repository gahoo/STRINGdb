test_get_homologs <- function() {
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=400 ) 
  data(diff_exp_example1)
  example1_mapped = string_db$map( diff_exp_example1, "gene", removeUnmappedRows = TRUE )  
  homologs = string_db$get_homologs( subset(example1_mapped, gene=="TP53")$STRING_id, target_species_id=10090  )
  checkTrue(nrow(homologs) > 0)
}
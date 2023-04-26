

#xxxxxxxxxxxxxx
# Changing node numbers to node names -------------------------------------
#xxxxxxxxxxxxxx
# nodenr_char: contains node numbers as character string
# tree_table: is a phylogenetic tree converted to tibble
func_nodenrs_to_nodenames <- function(nodenr_char, tree_table){
  nodenr_char %>% 
    as.numeric() %>% 
    tibble(node = .) %>% 
    left_join(., tree_table, by = "node") %>% 
    dplyr::select(label) %>% 
    pull()
}
# func_nodenrs_to_nodenames(nodenr_char = c(1, 2), tree_tab = tree_tab)
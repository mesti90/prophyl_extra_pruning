rm(list = ls())

library(tidyverse)
library(magrittr) # %<>% pipe
library(treeio) # as_tibble
library(phytools) # midpoint.root
library(ggtree) # tree plotting
library(Biostrings) # readDNAStringSet at func_alignment_length function
library(caper) # clade.members.list
library(knitr) # include_graphics
library(lubridate) # decimal_date

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- args[1]
  tree_path <- args[2]
} else {
  project_dir <- "~/Methods/prophyl"
  test_dir <- "~/Methods/prophyl-tests/test-shrink_tree"
  tree_path <- paste0(test_dir, "/chromosomes.node_labelled.final_tree.tre")
}

library(devtools)
load_all(project_dir)

# read trees calculated by Nextflow
tree <- ape::read.tree(tree_path)

# MANUAL: drop a tip
tree <- ape::drop.tip(tree, tip = "GCF_002992625.1")

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

# midpoint root
tree %<>% phytools::midpoint.root()

tree_tab <- tree %>% 
  # convert to tibble
  treeio::as_tibble() %>% 
  # adding a column: is.tip
  dplyr::mutate(is.tip = case_when(label %in% tree$tip.label ~ "Yes", TRUE ~ "No"))

# flag branches that are extremely long
# exclude zero branch length from calculating threshold
branch_length_qc_threshold <- qc(
  tree_tab$branch.length[which(tree_tab$branch.length > 1)],
  qlow = 0.25,
  qhigh = 0.75,
  m = 3,
  keep = "lower"
)$threshold[2]

tree_tab$branch.length.qc.pass <- tree_tab$branch.length < branch_length_qc_threshold 

# list of tip names of all internal nodes (defined with nrs.)
clade_members_list <- clade.members.list(
  phy = tree, 
  tip.labels = TRUE,
  include.nodes = TRUE) %>% 
  # converting a list of lists to a simple list
  list_flatten()

tree_tab$descendant_tips <- sapply(tree_tab$node, function(x) {
  index <- which(tree_tab$node == x)
  # if node is a tip, the number of descendant tips is 1.
  if (tree_tab$is.tip[index] == "Yes") {
    return(1) 
  # else the number of descendant tips is counted from clade_members_list
  } else {
    index2 <- which(names(clade_members_list) == paste0(x, "_tips"))
    if (length(index) == 0) {
      msg <- paste0(
        "Failed to count descendant tips for node: ",
        x
      )
      stop(msg)
    }
    length(clade_members_list[[index2]])
  }
})

#xxxxxxxxxx
# divide the list to 2 lists: tips & nodes
#xxxxxxxxxx
# Tips
#xxxxx
# clade_members_list_tips <- clade_members_list[grep(x = names(clade_members_list), 
#                                                    pattern = "_tips", 
#                                                    value = TRUE)]
# rename names: node nrs to labels
# names(clade_members_list_tips) %<>% 
#   # delete the "_tips" part from the name
#   gsub(x = ., pattern = "_tips", replacement = "") %>% 
#   as.numeric() %>% 
#   tibble(node = .) %>% 
#   left_join(., tree_tab, by = "node") %>% 
#   dplyr::select(label) %>% 
#   pull()

#xxxxx
# Nodes (the list contains node nrs. and not node names)
#xxxxx
# clade_members_list_nodes <- clade_members_list[grep(x = names(clade_members_list), 
#                                                    pattern = "_tips", 
#                                                    value = TRUE,
#                                                    invert = TRUE)]
# rm(clade_members_list)

# rename names: node nrs to labels
# names(clade_members_list_nodes) %<>%
#   # delete the "_nodes" part from the name
#   gsub(x = ., pattern = "_nodes", replacement = "") %>% 
#   as.numeric() %>% 
#   tibble(node = .) %>% 
#   left_join(., tree_tab, by = "node") %>% 
#   dplyr::select(label) %>% 
#   pull()

# change node nrs. to names in the list elements as well
# clade_members_list_nodes %<>% 
#   lapply(., func_nodenrs_to_nodenames, tree_tab = tree_tab)

# delete parts containing only a single node (the node itself as the list names)
# example: 
# $N263
# [1] "N263"
# -> delete elements which length is 1
# Nr_nodes <- sapply(clade_members_list_nodes, length)
# clade_members_list_nodes <- clade_members_list_nodes[names(Nr_nodes[Nr_nodes > 1])]
# rm(Nr_nodes)

# delete node names that are the same as the list names
# example:
# $N261
# [1] "N261" "N262" -> delete "N261" from the vector
# -> delete the first element in each vector
# clade_members_list_nodes %<>% 
#   lapply(., function(x) x[-1])

#xxxxxxxx
# nr of tips
#xxxxxxxx
# Nr_tips <- sapply(clade_members_list_tips, length)

# Merging tip nrs. to tree_tab_rooted
# tree_tab %<>% 
#   # add tip nrs.
#   left_join(., enframe(Nr_tips, name = "label", value = "Nr. of tips"), 
#             by = "label") %>% 
#   # add 1 as tip nr. when it is a tip
#   mutate(`Nr. of tips` = case_when(is.na(`Nr. of tips`) ~ 1,
#                                    TRUE ~ `Nr. of tips`))

# Inter Quartile Range: upper limit Q3 + 3 * (Q3 - Q1)
# qtr_vec <- tree_tab %>% 
#   filter(!is.na(branch.length)) %>% 
#   dplyr::select(branch.length) %>% 
#   pull() %>% 
#   quantile()

# iqr_bl_threshold <- qtr_vec["75%"] + 3 * (qtr_vec["75%"] - qtr_vec["25%"])
# names(iqr_bl_threshold) <- NULL
iqr_bl_threshold <- branch_length_qc_threshold

# Nr. of tips < %5 of all tips
nr_tips_threshold <- length(tree$tip.label)*0.05

# creating a tree_tab for plotting
tree_tab_plotting <- tree_tab %>% 
  # if branch.length is NA -> 0
  mutate(branch.length = case_when(is.na(branch.length) ~ 0,
                                   TRUE ~ branch.length))
# color by density
tree_tab_plotting$Density <- get_density(tree_tab_plotting$branch.length, 
                                         tree_tab_plotting$descendant_tips, 
                                         n = 100)
# scatterplot
g <- tree_tab_plotting %>% 
  ggplot(aes(x = branch.length, y = descendant_tips,
             color = Density)) +
  geom_point() +
  viridis::scale_color_viridis(option = "turbo") +
  xlab("Branch length") +
  # add branchlength threshold
  geom_vline(aes(xintercept = iqr_bl_threshold, 
                 linetype = "Branch length outliers\nto the right"), 
                 color = "darkred") +
  # add nr. of tips threshold
  geom_hline(aes(yintercept = nr_tips_threshold, 
                 linetype = "Nodes with less than 5%\nof the tips are below"),
                 color = "darkgreen") +
  # legend for lines
  scale_linetype_manual(name = "Outlier detection", values = c(2, 2), 
                        guide = guide_legend(override.aes = 
                                               list(color = c("darkred", 
                                                              "darkgreen")))) +
  theme_minimal()

if (!interactive()) {
  ggsave(
    filename = "outliers.pdf",
    plot = g,
    units = "cm",
    height = 8,
    width = 12
  )
}

# exclude tips that have extremely long terminal branches
exclude_tips <- tree_tab$label[which(
  tree_tab$is.tip == "Yes" &
  tree_tab$branch.length.qc.pass == FALSE
)]

# identify internal nodes that have extremely long brances
exclude_descendants <- tree_tab$node[which(
  tree_tab$is.tip == "No" &
  tree_tab$branch.length.qc.pass == FALSE &
  tree_tab$descendant_tips < nr_tips_threshold
)]

# exclude tips where an ancestral branch is extremely long
descendant_tips <- lapply(exclude_descendants, function(x) {
  index <- which(names(clade_members_list) == paste0(x, "_tips"))
  clade_members_list[[index]]
})

# highlight internal nodes where an ancestral branch is extremely long
# used for plotting
descendant_nodes <- lapply(exclude_descendants, function(x) {
  index <- which(names(clade_members_list) == paste0(x, "_nodes"))
  clade_members_list[[index]]
}) %>% unlist() %>% unique()

# concatenate the two vectors:
# tips that have extremely long terminal branches
# tips where an ancestral node has an extremely long branch
exclude_tips_all <- c(
  exclude_tips,
  unlist(descendant_tips)
) %>% unique()

tree_tab$exclude = FALSE
for (i in 1:nrow(tree_tab)) {
  if (tree_tab$label[i] %in% exclude_tips_all) {
    tree_tab$exclude[i] <- TRUE
  }
  if (tree_tab$node[i] %in% descendant_nodes) {
    tree_tab$exclude[i] <- TRUE
  }
}

# excluding those branches where length > iqr_bl_threshold & Nr. of tips < %5 of all tips

# vector of tips and internal nodes to exclude
# exc <- tree_tab %>% 
#   filter((branch.length > iqr_bl_threshold) & (`Nr. of tips` < nr_tips_threshold)) %>% 
#   dplyr::select(label) %>% 
#   pull()

# get the tips and further nodes based on exc
# exc <- clade_members_list_nodes[exc] %>% 
#   unlist(use.names = FALSE) %>% 
#   c(exc, .) %>% 
#   unique()

# exc <- clade_members_list_tips[exc] %>% 
#   unlist(use.names = FALSE) %>% 
#   c(exc, .) %>% 
#   unique()

# drop these tips to create a shrunken tree
shrinked_tree <- drop.tip(phy = tree, tip = exclude_tips_all)

# export shrinked tree

extension <- strsplit(tree_path, "\\.")[[1]]
extension <- extension[length(extension)]

if (!interactive()) {
  ape::write.tree(
    shrinked_tree,
    file = paste0("shrinked_tree.", estension)
  )
}

# tree plot
g2 <- tree_tab %>% 
  as.treedata() %>% 
  ggtree(aes(color = exclude), 
         layout = "equal_angle",
         size = 0.5) +
  scale_color_manual(values = c("black", "firebrick1")) +
  theme(legend.position = "bottom")

if (!interactive()) {
  ggsave(
    filename = "original_tree_marked_outliers.pdf",
    plot = g2,
    units = "cm",
    height = 20,
    width = 20
  )
}

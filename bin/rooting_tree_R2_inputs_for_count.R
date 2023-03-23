#!/usr/bin/env Rscript

# how to run from command line
# Rscript rooting_tree_R2_inputs_for_count_v1.2.R treeshrink.tre Coli_samples
# collection_date.tsv orthogroup_genome_prevalence_table.blastx.20230201.tsv.gz

#xxxxxxxxxxxxxxxx
# Installing required packages if needed, call them -----
#xxxxxxxxxxxxxxxx
start_time <- Sys.time()

library(tidyverse)
library(tidytree)
library(ape)
library(lubridate)
library(BactDating)
library(data.table)

#xxxxxxxxxxxxxxxxxxx
# Nextflow pipeline: after TreeShrink - before Count ------------------------
# Rooting with dates:
#   naming nodes
#   reading dates
#   calculating dated tree -> defining the root node
# creating root-2-tip plot and calculating R2 (with BactDating and with own script),
# creating newick tree file and prevalence tab tsv to Count software
#xxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxx
# Initialization ----------------------------------------------------------
#xxxxxxxxxxxxxxxxx
# Delete the files
outfile_names <- c("tree_rooting_R2_count_input.log",
                   "tree_rooted_by_time.tre",
                   "tree_dated.tre",
                   "bactdate_trace.pdf",
                   "bactdate_root_to_tip_regression.pdf",
                   "root_to_tip_result_tab.tsv",
                   "custom_root_to_tip_regression.pdf",
                   "for_count.tre",
                   "for_count_prevalence_tab.tsv")

unlink(outfile_names)
rm(outfile_names)

# creating log file
log_file <- "tree_rooting_R2_count_input.log"
# logging start time
sink(file = log_file, append = TRUE)
  cat(paste0("Start at: ", start_time, "\n"))
sink()

rm(start_time)

#xxxxxxxxxxxxxxxx
# * Checking the nr. of arguments ----
#xxxxxxxxxxxxxxxx
args = commandArgs(trailingOnly = TRUE)

# test if there is at least 3 arguments: if not, return an error
w <- paste0(
  "ERROR: Please supply the following files as arguments in the following ",
  "order: a newick tree, a tsv with dates, a gene prevalence table in tsv or ",
  "tsv.gz, an SNP alignment.\n"
)
if (length(args) != 4) {
  sink(file = log_file, append = TRUE)
  cat(w)
  sink()
  stop(w, call. = FALSE)
}
rm(w)

#xxxxxxxxxxxxxxxx
# Tree --------------------------------------------------------------------
#xxxxxxxxxxxxxxxx
# read tree file
tree <- read.tree(args[1]) %>% 
  # transforms all multichotomies into a series of dichotomies
  multi2di() %>% 
  # unroot tree
  unroot()

#xxxxxxxxxx
# * Add node names to tree ----
#xxxxxxxxxx
# convert tree to tibble
tree_tab <- tree %>% 
  as_tibble() %>% 
  # create Node names
  mutate(label = case_when(is.na(label) ~ paste0("N", node),
                           TRUE ~ label))

# convert branch length from per site to per genome
f <- seqinr::read.fasta(args[4])
alignment_length <- length(f[[1]])
tree_tab$branch.length <- tree_tab$branch.length*alignment_length

# convert tibble back to phy
tree <- tree_tab %>% 
  as.phylo()

rm(tree_tab)

#xxxxxxxxxxxxxxxxxxx
# Date the tips -------------------------------------------
#xxxxxxxxxxxxxxxxxxx

# import dates for tips
date_tab <- read_tsv(args[2], col_names = FALSE) %>% 
  # rename columns
  dplyr::rename(label = X1, `Sample collection date` = X2) %>% 
  # keep tips of the actual tree only
  filter(label %in% tree$tip.label & !is.na(`Sample collection date`)) %>% 
  # create a new column with decimal dates
  mutate(`Sample collection date decimal` = decimal_date(`Sample collection date`))

# add tip labels with missing dates
date_tab <- tree$tip.label %>% 
  as_tibble() %>% 
  rename(label = value) %>% 
  # arranged by coli_tree$tip.label
  full_join(., date_tab)

#xxxxxxxxxxxxxxxxxxxx
# Rooting with dated tips ------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxx

# standard argument settings
# (tree, date, initMu = NA, initAlpha = NA, initSigma = NA, 
# updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T, 
# nbIts = 10000 Nr. of MCMC iterations
# thin = ceiling(nbIts/1000), useCoalPrior = T, 
# model = "arc", useRec = F, minbralen = 0.1, showProgress = F) 

dated_tree <- bactdate(tree = tree, 
                       date = date_tab$`Sample collection date decimal`,
                       model = "arc")

# now, the dated_tree$inputtree contains the root
# overwrite the original tree
tree <- dated_tree$inputtree

# write trees to files
tree %>% write.tree("tree_rooted_by_time.tre")
dated_tree$tree %>% write.tree("tree_dated.tre")

#xxxxxxx
# * BactDate MCMC trace plot ----
#xxxxxxx
pdf("bactdate_trace.pdf", width = 10, height = 6)
  plot(dated_tree, 'trace')
dev.off()

#xxxxxxxx
# * BactDate root-to-tip plot ----
#xxxxxxxx
pdf("bactdate_root_to_tip_regression.pdf", width = 10, height = 6)
  roottotip(tree, date_tab$`Sample collection date decimal`)
dev.off()

#xxxxxxxxx
# Root-to-tip distance ----------------------------------------------------
#xxxxxxxxx
root_to_tip_result_tab <- node.depth.edgelength(tree) %>% 
  as_tibble() %>% 
  rename("root-to-tip-distance" = "value") %>% 
  # add tip and node names
  bind_cols(label = c(tree$tip.label, tree$node.label), .) %>% 
  # add collection date
  left_join(., date_tab, by = "label") %>% 
  # add a column containing if to plot
  mutate(to_plot = !is.na(`Sample collection date decimal`))

# write Root-to-tip distances to tsv
root_to_tip_result_tab %>% 
  write_tsv("root_to_tip_result_tab.tsv")

#xxxxxxxxxx
# * Calculating R2 of the linear regression ----
#xxxxxxxxxx
Lin_reg_result <- root_to_tip_result_tab %>% 
  filter(to_plot == TRUE) %>% 
  lm(`root-to-tip-distance` ~ `Sample collection date decimal`, data = .) %>% 
  summary()

sink(file = log_file, append = TRUE)
  cat(paste0("DATA: Adjusted R2 of linear regression: ", Lin_reg_result$adj.r.squared, "\n"))
sink()

#xxxxxxxxxxxxx
# Plot root to tip regression --------------------------------------------------
#xxxxxxxxxxxxx
p <- root_to_tip_result_tab %>% 
  filter(to_plot == TRUE) %>% 
  ggplot(aes(x = `Sample collection date decimal`,
             y = `root-to-tip-distance`)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "steelblue") +
  # name outliers
  # geom_text_repel(aes(label = label), size = 3, max.overlaps = 25) +
  # rename x axis
  xlab("Sample collection date") +
  # add title with R2
  ggtitle(paste0("Adjusted R2: ", round(Lin_reg_result$adj.r.squared, 4))) +
  theme_bw()

ggsave(plot = p, filename = "custom_root_to_tip_regression.pdf", 
       width = 8, height = 8)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Preparing files for count -----------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxx
# * Tree: enshort tip names for count ----
#xxxxxxxxxx
# convert tree to tibble
tree_tab <- tree %>% 
  as_tibble() %>% 
  # convert "Coli_" to "C" in tip names
  mutate(label = case_when(!is.na(label) ~ str_remove(string = label, 
                                                      pattern = "oli_"),
                           TRUE ~ label))

# convert tibble to phy
tree <- tree_tab %>% 
  as.phylo()

tree %>% write.tree("for_count.tre")

#xxxxxxxxxx
# * Orthogroup prevalence table ----
#xxxxxxxxxx
fread(args[3]) %>% 
  select(c("Custom_ID", all_of(tree$tip.label))) %>% 
  write_tsv("for_count_prevalence_tab.tsv")

# logging end time
sink(file = log_file, append = TRUE)
  cat(paste0("The tree rooting, count input file preparing code has finished at: ", Sys.time(), "\nSession info:\n"))
  sessionInfo()
sink()

quit(save = "yes")

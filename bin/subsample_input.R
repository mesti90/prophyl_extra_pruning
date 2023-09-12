# This script takes a table of assemblies and draws samples using the selected
# sampling strategy. Each subsampled table will then be used to build and date a
# subset tree. Each subset tree will be used to calculate a number of data
# points for the relative risk plots.

library(optparse)
rm(list=ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("A", "--ASSEMBLIES"),
    type = "character",
    help = "Path to the assemblies file that was used for building the tree."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to assemblies file used for subsampling."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
  ),
  make_option(
    c("-d", "--duplicates"),
    type = "character",
    help = "A text file which contains tip labels for identical tips."
  ),
  make_option(
    c("-c", "--subsample_count"),
    type = "character",
    help = "Number of subsample sets to draw."
  ),
  make_option(
    c("-C", "--subsample_tipcount"),
    type = "character",
    help = "Number of tips to draw in each subsample set."
  ),
  make_option(
    c("-b", "--focus_by"),
    type = "character",
    help = "Variable to focus on."
  ),
  make_option(
    c("-o", "--focus_on"),
    type = "character",
    help = "Value of the variable to focus on."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    ASSEMBLIES = "assemblies.tsv",
    assemblies = "assemblies_for_country_rr.tsv",
    tree = "dated_tree.rds",
    duplicates = "duplicates.txt",
    subsample_count = 10,
    subsample_tipcount = 250,
    focus_by = "continent",
    focus_on = "europe"
  )
}

library(devtools)
load_all(args$project_dir)

# read assemblies used for building tree - required for managing duplicates
ASSEMBLIES <- read.csv(args$ASSEMBLIES, sep = "\t")
# read assemblies - this is the table which will be subsampled
assemblies <- read.csv(args$assemblies, sep = "\t")
# read tree
tree <- readRDS(args$tree)
# import duplicates
duplicates <- parse_duplicates(args$duplicates)
# number of subsample sets to draw
subsample_count <- as.numeric(args$subsample_count)
# number of tips to draw in each subsample set
subsample_tipcount <- as.numeric(args$subsample_tipcount)

# import focus variables
# a variable within the input table used for focusing.
# only used if sampling strategy is "focused".
focus_by <- args$focus_by
focus_on <- args$focus_on

if (focus_by == "none") focus_by <- NULL
if (focus_on == "none") focus_on <- NULL

focus <- validate_focus(assemblies, focus_by, focus_on)
focus_by <- focus$focus_by
focus_on <- focus$focus_on

# ratio of samples to take from focus group.
focus_ratio <- 0.75

# define sampling strategy, "random", or "focused"
if (is.null(focus_by)) {
  type <- "random"
} else {
  type <- "focused"
}

# The shrinked tree may contain less tips than the original tree
# Only sample assemblies that are included in the shrinked tree
index <- which(assemblies$assembly %in% tree$tip.label == FALSE)
if (length(index) > 0) {
  assemblies <- assemblies[-index, ]
}

# The number of tips in each subsample must be smaller than the overall number
# of tips. If not, stop with an informative error.
if (subsample_tipcount >=  length(tree$tip.label)) {
  stop(paste0(
    "Parameter 'subsample_tipcount' must be lower than ",
    "the overall number of tips."
  ))
}

subs <- list()
if (type == "random") {
  set.seed(0)
  for (i in 1:subsample_count) {
    subs[[i]] <- assemblies[sample(1:nrow(assemblies), subsample_tipcount, replace = FALSE), ]
  }
}

if (type == "focused") {
  set.seed(0)
  focus_count <- round(focus_ratio * subsample_tipcount, 0)
  no_focus_count <- subsample_tipcount - focus_count
  if (focus_count == 0 | no_focus_count == 0) {
    stop("Number of assemblies in and outside of focus group must be non zero.")
  }
  for (i in 1:subsample_count) {
    focus_index <- which(assemblies[[focus_by]] == focus_on)
    if (length(focus_index) >= focus_count) {
      focus <- sample(
      focus_index,
      focus_count,
      replace = FALSE
      )
    } else {
      stop("Not enough assemblies in focus group to subsample")
    }
    # exclude any assemblies where focus variable is NA
    no_focus_index <- which(assemblies[[focus_by]] != focus_on)
    if (length(no_focus_index) >= no_focus_count) {
      no_focus <- sample(
      no_focus_index,
      no_focus_count,
      replace = FALSE
      )
    } else {
      stop("Not enough assemblies in non-focus group to subsample")
    }
    subs[[i]] <- assemblies[c(focus, no_focus), ]
  }
}

# manage duplicates

subs_tidy <- list()
dups <- list()
for (i in 1:subsample_count) {
  tidydbs <- tidy_duplicates(
    df = ASSEMBLIES,
    subset = subs[[i]], 
    id_var = "assembly",
    duplicates = duplicates
  )
  
  subs_tidy[[i]] <- tidydbs$subset
  dups[[i]] <- tidydbs$duplist
  
  subset <- tidydbs$subset
  duplist <- tidydbs$duplist
  
  # consistency checks
  
  # all assemblies have been replaced with reference assemblies where necessary
  # this is relevant when there are duplicates in the assembly data set
  if (length(duplicates) > 0 & length(duplist) > 0) {
    dupnames <- unname(unlist(duplist))
    refnames <- names(duplist)
    # these are the refs that are also in subset
    refs_in_sub <- refnames[which(refnames %in% dupnames)]
    # these are the duplicates that were removed
    dups_removed <- dupnames[-which(dupnames %in% refs_in_sub)]
    # these are the refs that were added
    refs_added <- refnames[which(!refnames %in% refs_in_sub)]
    # refs are really in subset
    testthat::expect_true(all(refs_in_sub %in% subset$assembly))
    testthat::expect_true(all(refs_added %in% subset$assembly))
    # removed dups are not in subset
    testthat::expect_false(any(dups_removed %in% subset$assembly))
    # subset has correct number of rows
    removed_count <- sum(sapply(duplist, function(x) length(x)-1))
    testthat::expect_equal(
      nrow(subset),
      nrow(subs[[i]])-removed_count
    )
    # subset still has unique rows
    testthat::expect_equal(nrow(subset), length(unique(subset$assembly)))
  }
}

for (i in 1:subsample_count) {
  digits <- ceiling(log10(subsample_count+1))
  zeroes <- digits - floor(log10(i))
  filename = paste0(
    c("subsample_", rep(0, times = zeroes), i, ".tsv"), collapse = "")
  dupfile = paste0(
    c("subsample_", rep(0, times = zeroes), i, ".rds"), collapse = "")
  # export subset
  write.table(
    subs_tidy[[i]],
    file = filename,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  # export duplicates
  saveRDS(dups[[i]], dupfile)
}

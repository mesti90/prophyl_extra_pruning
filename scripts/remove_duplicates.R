# TODO
# This script is very inefficient because it requires me to read in all the 
# chromosomes. This can easily lead to memory issues. Therefore it is better
# to use less memory intensive tools to just remove duplicates and print a 
# list of genomes that were removed and later rename some tips if necessary.

library(dplyr)
library(optparse)
rm(list = ls())

con <- file("log.txt")
sink(con, split = TRUE)

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to the assemblies file."
  ),
  make_option(
    c("-f", "--fasta_file"),
    type = "character",
    help = "Path to the fasta file containing pseudo-whole chromosome sequences."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    assemblies = "aci_all.tsv",
    fasta_file = "chromosomes.fasta"
  )
}

# read fasta file with pseudo-whole chromosome sequences
f <- Biostrings::readDNAStringSet(args$fasta_file)

# this does not track first occurrence but tracks subsequent occurrences
index_duplicated <- which(duplicated(f))

# read assemblies, required for finding earliest sample within multiplicates
assemblies <- read.csv(args$assemblies, sep = "\t")

# create a list where each element is a vector of assembly names that have identical sequence
# the name of each element should be the name of the assembly which will be kept.
multiplicates <- list()
if (length(index_duplicated) > 0) {
  for (i in seq_along(index_duplicated)) {
    index <- which(f == f[index_duplicated][i])

    multiplicates[[i]] <- names(f)[index]

    # get collection dates from metadata
    dates <- sapply(multiplicates[[i]], function(x) {
        index <- which(assemblies$assembly == x)
        if (length(index) == 0) {
          warning(paste0("Not found: ", x))
          return(NA)
        }
        if (length(index) == 1) {
          return(assemblies$collection_date[index])
        }
        if (length(index) > 1) {
          stop(paste0("Found multiple times in database: ", x))
        }
    })

    # if all dates are NA, keep the first in the list
    if (all(is.na(dates))) {
        names(multiplicates)[i] <- multiplicates[[i]][1]
        next()
    }

    # find earliest sample by year
    years <- sapply(multiplicates[[i]], function(x) {
        index <- which(assemblies$assembly == x)
        if (length(index) == 0) {
          warning(paste0("Not found: ", x))
          return(NA)
        }
        if (length(index) == 1) {
          return(assemblies$collection_year[index])
        }
        if (length(index) > 1) {
          stop(paste0("Found multiple times in database: ", x))
        }
    }) %>% as.numeric()

    # if all years are NA, keep the first in the list
    if (all(is.na(years))) {
        names(multiplicates)[i] <- multiplicates[[i]][1]
        next()
    }

    # if the earliest sample is a single sample, great
    index_year <- which(years == min(years, na.rm = TRUE))
    if (length(index_year) == 1) {
        names(multiplicates) <- names(f)[index][index_year]
        next()
    }

    # if not, find the earliest sample by day
    # TODO wrap this into try and return NA if cannot be converted
    days <- as.Date(dates[index_year], origin = "1970-01-01")

    # if all day are NA, keep the first in the list
    if (all(is.na(days))) {
        names(multipcates)[i] <- names(f)[index][index_year][1]
        next()
    }

    # if the earliest sample is a single sample, great
    # if there are more than one, keep the first
    index_day <- which(day == min(days, na.rm = TRUE))
    names(multiplicates)[i] <- names(f)[index][index_year][index_day][1]    
  }
  
  for (i in 1:length(multiplicates)) {
    multiplicates[[i]] <- data.frame(
      include = names(multiplicates)[i],
      exclude = multiplicates[[i]]
    )
  }
  
  multiplicates <- dplyr::bind_rows(multiplicates)
  
  # remove rows where "include" and "exclude" are the same
  multiplicates <- multiplicates[which(multiplicates$include != multiplicates$exclude),]
  
  multiplicates$include_date <- sapply(multiplicates$include, function(x) {
    index <- which(assemblies$assembly == x)
    if (length(index) == 0) {
      warning(paste0("Not found: ", x))
      return(NA)
    }
    if (length(index) == 1) {
      return(assemblies$collection_date[index])
    }
    if (length(index) > 1) {
      stop(paste0("Found multiple times in database: ", x))
    }
  })
  multiplicates$exclude_date <- sapply(multiplicates$exclude, function(x) {
    index <- which(assemblies$assembly == x)
    if (length(index) == 0) {
      warning(paste0("Not found: ", x))
      return(NA)
    }
    if (length(index) == 1) {
      return(assemblies$collection_date[index])
    }
    if (length(index) > 1) {
      stop(paste0("Found multiple times in database: ", x))
    }
  })
  
  # check that none of the "include" entries equal an "exclude" entry
  testthat::expect_true(all(multiplicates$include != multiplicates$exclude))
  
  # remove sequences that were marked from removal
  index_remove <- unique(which(names(f) %in% multiplicates$exclude))
  if (length(index) > 0) {
    f_unique <- f[-index_remove]
  }
  
  # check that sequences are no longer duplicated
  testthat::expect_equal(sum(duplicated(f_unique)), 0)
  
  # check that sequences which should not be excluded are still there
  testthat::expect_true(all(multiplicates$include %in% names(f_unique)))
  
  # check that sequences which should be excluded are no longer there
  testthat::expect_true(all(!multiplicates$exclude %in% names(f_unique)))
}

if (length(index_duplicated) > 0) {
  # export deduplicated fasta
  Biostrings::writeXStringSet(f_unique, filepath = "chromosomes_dedup.fasta")
} else {
  cmd = paste0("ln -s ", args$fasta_file, " chromosomes_dedup.fasta")
  system(cmd)
}

# export list of sequences that were removed
saveRDS(multiplicates, file = "multiplicates.rds")
  
# close log file
sink(con)
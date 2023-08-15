library(dplyr)
library(ggplot2)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-t", "--trees"),
    type = "character",
    help = "All dated trees collected in a single rds file."
  ),
  make_option(
    c("-s", "--sensitivity"),
    type = "character",
    help = "A table of phage-host sensitivity test results."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    trees = "all_dated_trees.rds",
    sensitivity = "Acinetobacter baumannii strains_Phages - Acinetobacter baumannii strains_Phages.tsv"
  )
}

# import trees
trees <- readRDS(args$trees)

# calculate cophenetic (patristic) distance matrix for each pair of genomes
patrdist <- lapply(trees, ape::cophenetic.phylo)
names(patrdist) <- names(trees)

# import sensitivity test results
phres <- read.csv(args$sensitivity, sep = "\t", na.strings = c("", "NA"))

# only keep strains where sequences are available
phres <- phres[which(!is.na(phres$Assembly)),]

# only keep strains where tree for the MLST is available
phres <- phres[which(phres$MLST %in% names(trees)),]

# only keep strains which were tested against all phages
phres <- phres[complete.cases(phres[,8:ncol(phres)]),]

# only keep strains which are present on a tree
# only keep strains which are present on a single tree
remove <- vector()
for (i in phres$Assembly) {
  hit <- FALSE
  for (j in names(trees)) {
    if (i %in% trees[[j]]$tip.label) {
      if (hit == FALSE) {
        hit <- TRUE
      } else {
        msg <- paste0(i, " found on multiple trees!")
        message(msg)
      }
    }
  }
  if (hit == FALSE) {
    msg <- paste0(i, " not found!")
    message(msg)
    remove <- c(remove, i)
  }
}

if (length(remove) > 0) {
  phres <- phres[-which(phres$Assembly %in% remove),]
}

# format phres
for (i in 1:nrow(phres)) {
  for (j in 8:ncol(phres)) {
    # Count partial sensitivity to a phage as no sensitivity
    if (!is.na(phres[i,j]) && phres[i,j] == "1") {
      phres[i,j] <- "0"
    }
    if (!is.na(phres[i,j]) && phres[i,j] %in% c("YES", "Yes")) {
      phres[i,j] <- "1"
    }
    if (!is.na(phres[i,j]) && phres[i,j] == "NO") {
      phres[i,j] <- "0"
    }
  }
}

for (j in 8:ncol(phres)) {
  phres[,j] <- as.numeric(phres[,j])
  phres[,j] <- ifelse(
    phres[,j] > 1,
    1,
    phres[,j]
  )
}

# only keep strains which were sensitive to at least one phage
# effective_phage_count <- apply(phres[,8:ncol(phres)], 1, sum)
# phres <- phres[which(effective_phage_count > 0),]

# only keep strains that belong to KL types against which at least one phage
# was found. e.g. if there was an effective phage against KL3 then keep all KL3
# strains, regardless of their phage profile.

infection_count <- sapply(unique(phres$KL), function(x) {
  sum(phres[which(phres$KL == x),8:ncol(phres)])
})
kl_to_keep <- names(infection_count)[which(infection_count > 0)]
phres <- phres[which(phres$KL %in% kl_to_keep),]

comp <- combn(phres$Assembly, 2) %>% t() %>% as.data.frame()
names(comp) <- c("A1", "A2")
comp$A1_mlst <- sapply(comp$A1, function(x) {
  index <- which(phres$Assembly == x)
  phres$MLST[index]
})
comp$A2_mlst <- sapply(comp$A2, function(x) {
  index <- which(phres$Assembly == x)
  phres$MLST[index]
})
comp$A1_KL <- sapply(comp$A1, function(x) {
  index <- which(phres$Assembly == x)
  phres$KL[index]
})
comp$A2_KL <- sapply(comp$A2, function(x) {
  index <- which(phres$Assembly == x)
  phres$KL[index]
})
# only compare strains which had the same MLST and K type
comp <- comp[which(comp$A1_mlst == comp$A2_mlst & comp$A1_KL == comp$A2_KL),]

# copy patristic distances from distance matrix
comp$patristic <- NA
for (i in 1:nrow(comp)) {
  ST <- comp$A1_mlst[i]
  index_x <- which(rownames(patrdist[[ST]]) == comp$A1[i])
  index_y <- which(colnames(patrdist[[ST]]) == comp$A2[i])
  comp$patristic[i] <- patrdist[[ST]][index_x, index_y]
}

# since distance is in years, divide by two to relate to divergence times.
comp$patristic <- comp$patristic/2

# calculate phage sensitivity distance matrix
phagedist <- as.matrix(vegan::vegdist(phres[8:ncol(phres)], "jaccard", na.rm = FALSE))
rownames(phagedist) <- phres$Assembly
colnames(phagedist) <- phres$Assembly

# copy phage sensitivity distances from distance matrix
comp$phagedist <- NA
for (i in 1:nrow(comp)) {
  index_x <- which(rownames(phagedist) == comp$A1[i])
  index_y <- which(colnames(phagedist) == comp$A2[i])
  comp$phagedist[i] <- phagedist[index_x, index_y]
}

comp$patristic <- signif(comp$patristic, 4)
comp$phagedist <- signif(comp$phagedist, 4)

write.table(
  comp,
  file = "phylodist_phagedist.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

comp_jitter <- comp

set.seed(0)
comp_jitter$patristic <- jitter(comp$patristic)
comp_jitter$phagedist <- jitter(comp$phagedist)

# define combined serotype to use with facet wrap
comp_jitter$serotype <- paste0(comp_jitter$A1_mlst, " - ", comp_jitter$A1_KL)

# filter to serotypes with at least 15 comparisons (6 strains)
serotypes <- comp_jitter %>% 
  group_by(serotype) %>% 
  summarise(count = n()) %>% 
  filter(count >= 15)
comp_jitter <- comp_jitter[which(comp_jitter$serotype %in% serotypes$serotype),]

# define limits for plotting
xmin <- min(comp_jitter$patristic, na.rm = TRUE)
xmax <- max(comp_jitter$patristic, na.rm = TRUE)
ymin <- min(comp_jitter$phagedist, na.rm = TRUE)
ymax <- max(comp_jitter$phagedist, na.rm = TRUE)

g <- ggplot(comp_jitter, aes(patristic, phagedist)) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", span = 5) +
  xlim(xmin, xmax) +
  ylim(ymin, ymax) +
  xlab("Time to most recent common ancestor (years)") + 
  ylab("Euclidean distance of phage sensitivity") +
  facet_wrap(serotype~.)

ggsave(
  filename = "phylodist_phagedist.pdf",
  plot = g,
  units = "cm",
  height = 20,
  width = 20
)

# ANALYSE PAIRS WHICH HAD IDENTICAL PHAGE SENSITIVITY PROFILES

# Assign strains with identical phage sensitivity profiles to the same cluster
d <- dist(phres[,8:ncol(phres)])
h <- hclust(d)
phres$cluster <- cutree(h, h = 0)

# tidy cluster names
nchar_max <- max(nchar(phres$cluster))
phres$cluster <- sapply(phres$cluster, function(x) {
  x <- as.character(x)
  n <- nchar(x)
  zeros <- nchar_max - n
  if(zeros == 0) {
    out2 <- paste0("cluster_", x)
  }
  if (zeros > 0) {
    out <- strsplit(x, "")[[1]]
    out2 <- rep(NA, times = nchar_max + 8)
    out2[1:8] <- c("c","l","u","s","t","e","r","_")
    out2[9:(9+zeros-1)] <- "0"
    out2[(9+zeros):length(out2)] <- out
    out2 <- paste(out2, collapse = "")
  }
  return(out2)
})

write.table(
  phres,
  file = "phage_host_resistance_profiles.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

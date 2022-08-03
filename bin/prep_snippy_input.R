rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

con <- file("log.txt", open = "wt")
sink(con, type = "message")

df <- read.csv(args[1], sep = "\t")
source_dir <- args[2]

paired_reads <- data.frame()
single_reads <- data.frame()
contigs <- data.frame()

for (i in 1:nrow(df)) {
  message(paste0("Searching ", df$assembly[i], ". "), appendLF = FALSE)
  read_dir <- paste0(source_dir, "/", df$jobname[i], "/raw_reads/")
  if (dir.exists(read_dir)) {
    message("Raw reads directory found. ", appendLF = FALSE)
    hit <- list.files(read_dir)[grep(df$assembly[i], list.files(read_dir))]
    if (length(hit) == 0) message(" Raw reads not found. Skipping.")
    if (length(hit) == 1) {
      message(" Single reads found.")
      new_row <- data.frame(
        assembly = df$assembly[i],
        reads = paste0(read_dir, hit)
        )
      single_reads <- rbind(single_reads, new_row)
    }
    if (length(hit) > 1) {
      d <- stringdist::stringdist(hit, df$assembly[i])
      index <- which(d == min(d))
      if (length(index) == 1) {
        message(" Single reads found.")
        new_row <- data.frame(
          assembly = df$assembly[i],
          reads = paste0(read_dir, hit[index])
        )
        single_reads <- rbind(single_reads, new_row)
      }
      if (length(index) == 2) {
        message(" Paired reads found.")
        index.R1 <- grep("R1", hit[index])
        index.R2 <- grep("R2", hit[index])
        new_row <- data.frame(
          assembly = df$assembly[i],
          R1 = paste0(read_dir, hit[index][index.R1]),
          R2 = paste0(read_dir, hit[index][index.R2])
        )
        paired_reads <- rbind(paired_reads, new_row)  
      }
      if (length(index) > 2) {
        message(paste0(" More than 2 files found. Skipping."))
      }
    }
  } else {
    assembly_dir <- paste0(source_dir, "/", df$jobname[i], "/assembled_genomes/")
    if (dir.exists(assembly_dir)) {
      message("Assembly directory found. ", appendLF = FALSE)
      hit <- list.files(assembly_dir)[grep(df$assembly[i], list.files(assembly_dir))]
      if (length(hit) == 0) message(" Assembly not found. Skipping.")
      if (length(hit) == 1) {
        message(" Assembly found.")
        new_row <- data.frame(
          assembly = df$assembly[i],
          path = paste0(assembly_dir, hit)
        )
        contigs <- rbind(contigs, new_row)
      }
      if (length(hit) > 1) {
        d <- stringdist::stringdist(hit, df$assembly[i])
        index <- which(d == min(d))
        if (length(index) == 1) {
          message(" Assembly found.")
          new_row <- data.frame(
            assembly = df$assembly[i],
            path = paste0(assembly_dir, hit[index])
          ) 
          contigs <- rbind(contigs, new_row)
        } else {
          message(paste0(" More than 1 files found. Skipping."))
        }
      }
    } else {
      message(paste0("Not found. Skipping."))
    }
  }
}

if (nrow(paired_reads) == 0) {
  paired_reads <-  paste("assembly", "R1", "R2", sep = "\t")
  write.table(
    paired_reads,
    file = "paired_reads.csv",
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE)
} else {
  write.csv(
    paired_reads, file = "paired_reads.csv", row.names = FALSE, quote = FALSE)
}

if (nrow(single_reads) == 0) {
  single_reads <-  paste("assembly", "reads", sep = "\t")
  write.table(
    single_reads,
    file = "single_reads.csv",
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE)
} else {
  write.csv(
    single_reads, file = "single_reads.csv", row.names = FALSE, quote = FALSE)
}

if (nrow(contigs) == 0) {
  contigs <- paste("assembly", "path", spe = "\t")
  write.table(
    contigs,
    file = "contigs.csv",
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE)
} else {
  write.csv(
    contigs, file = "contigs.csv", row.names = FALSE, quote = FALSE)
}

sink(con)


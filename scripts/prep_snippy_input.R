rm(list=ls())
con <- file(paste0(
  "./jobs/", jobname, "/output/","log_", format(Sys.time(), "%Y%m%d"), ".txt"))
sink(con, append = TRUE, split = TRUE)

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]
source_dir <- args[2]

output_dir <- paste0("./jobs/", jobname, "/output/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

df <- read.csv(paste0("./jobs/", jobname, "/assemblies.tsv"),
               sep = "\t")

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

write.table(paired_reads,
            file = paste0(output_dir, "paired_reads.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")


write.table(single_reads,
            file = paste0(output_dir, "single_reads.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(contigs,
            file = paste0(output_dir, "contigs.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

sink()
close.connection(con)

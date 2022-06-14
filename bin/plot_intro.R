#rm(list = ls())
library(dplyr)

# dates from Hun to Eng format
Sys.setlocale("LC_TIME", "C")

#args <- commandArgs(trailingOnly = TRUE)
#tree_tbl_file <- args[1]

#tree_tbl <- readRDS(tree_tbl_file)

# only keep rows where only one country was predicted for the destination node
index <- grep("\\|", tree_tbl$country)
if (length(index) > 0) {
  tree_tbl <- tree_tbl[-index,]
}

for (i in unique(tree_tbl$country)) {
  index <- which(tree_tbl$country == i)
  country_tbl <- tree_tbl[index, ]
  
  country_sum <- country_tbl %>%
    group_by(collection_year_pred) %>%
    summarise(count = sum(p_intro))
  
  g <- ggplot(country_sum, aes(collection_year_pred, count))+
    geom_col()+
    ylab("Number of introductions")+
    ggtitle(i)
  
  ggsave(filename = paste0(i, ".pdf"))
}  


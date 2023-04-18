FROM debian:testing-20230411

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# Install R and R packages from apt
RUN apt-get install -y \
  r-base \
  r-bioc-biostrings \
  r-cran-ape \
  r-cran-biocmanager \
  r-cran-data.table \
  r-cran-devtools \
  r-cran-dplyr \
  r-cran-ggrepel \
  r-cran-ggplot2 \
  r-cran-hmisc \
  r-cran-htmlwidgets \
  r-cran-igraph \
  r-cran-lubridate \
  r-cran-magrittr \
  r-cran-reshape2 \
  r-cran-rjava \
  r-cran-seqinr \
  r-cran-stringdist \
  r-cran-stringi \
  r-cran-tibble \
  r-cran-tidyr \
  r-cran-tidyverse

# Install ggmap
RUN R -e 'options(warn = 2); install.packages("cli", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("ggmap", dependencies = TRUE)'

# Install leaflet
RUN apt-get install -y \
  libgdal-dev \
  libudunits2-dev
RUN R -e 'options(warn = 2); install.packages("leaflet", dependencies = TRUE)'

# Install R packages from CRAN
RUN R -e 'options(warn = 2); install.packages("leaflet.minicharts", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("qualpalr", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("tictoc", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("tidytree", dependencies = TRUE)'

# Install procps 
# https://github.com/replikation/What_the_Phage/issues/89
RUN apt-get install -y procps

# Install geosphere
RUN apt-get install -y r-cran-geosphere

# Install treedater from GitHub instead of CRAN
RUN R -e 'options(warn = 2); devtools::install_github("emvolz/treedater@7b8a72aa0ea71ded4cde6cd287529aa8ae679c68")'

# Install BactDating
RUN R -e 'options(warn = 2); devtools::install_github("xavierdidelot/BactDating@v1.1")'

# Install phytools
RUN apt-get install -y r-cran-phytools

# Install ggtree
RUN R -e 'options(warn = 2); BiocManager::install("ggtree", update = FALSE, ask = FALSE)'

# Install ggimage, ggnewscale, ggtreeExtra
RUN R -e 'options(warn = 2); install.packages("ggimage", dependencies = TRUE)'
RUN R -e 'options(warn = 2); install.packages("ggnewscale", dependencies = TRUE)'
RUN R -e 'options(warn = 2); BiocManager::install("ggtreeExtra", update = FALSE, ask = FALSE)'

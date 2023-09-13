FROM debian:bookworm-20230904

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

# # Install leaflet
# RUN apt-get install -y \
#   libgdal-dev \
#   libudunits2-dev
# RUN R -e 'options(warn = 2); install.packages("leaflet", dependencies = TRUE)'

# Install R packages from CRAN
# RUN R -e 'options(warn = 2); install.packages("leaflet.minicharts", dependencies = TRUE)'
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

# # Install quarto
# RUN apt-get install -y wget
# RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.3.336/quarto-1.3.336-linux-amd64.deb
# RUN dpkg -i quarto-1.3.336-linux-amd64.deb

# Install optparse
RUN apt-get install -y r-cran-optparse

# Install Gubbins dependencies - FastTree, IQTree, RAxML
RUN apt-get install -y fasttree iqtree raxml 

RUN mkdir iqtree2 && ln -s /bin/iqtree2 /iqtree2/iqtree
ENV PATH="${PATH}:/iqtree2"

# Install Gubbins dependencies - RAxML-NG 
RUN apt-get install -y bison clang-tidy cmake flex git libgmp3-dev libgtest-dev make mpi-default-dev

RUN git clone --recursive https://github.com/amkozlov/raxml-ng && \
  cd raxml-ng && \
  mkdir build && \
  cd build && \
  cmake .. && \
  make

ENV PATH="${PATH}:/raxml-ng/bin"

# Install Gubbins dependencies - Rapidnj
RUN git clone https://github.com/somme89/rapidNJ.git && \
  cd rapidNJ && \
  make -j 8 && \
  mv bin/rapidnj /usr/bin/rapidnj && \
  cd .. && \
  rm -rf rapidNJ

# Install Gubbins dependencies - Python modules
RUN apt-get install -y python3 python3-dev python3-pip python3-tk python3-yaml

RUN apt-get install -y \
  python3-biopython \
  python3-dendropy \
  python3-multiprocess \
  python3-numba \
  python3-numpy \
  python3-scipy

# Install Gubbins dependencies
RUN apt-get install -y \
  autoconf-archive \
  check \
  curl \
  gcc \
  libmpich-dev \
  libtool \
  mpich

RUN git clone https://github.com/nickjcroucher/gubbins && \
  cd gubbins && \
  git checkout 21e70def99edad2d47c150c3ed291ad11350a2d3 && \
  autoreconf -i && \
  ./configure && make && make install

ENV OPENBLAS_NUM_THREADS=2
ENV LD_LIBRARY_PATH="/usr/local/lib/:$LD_LIBRARY_PATH"

# Install TreeTools
RUN R -e 'options(warn = 2); install.packages("TreeTools", dependencies = TRUE)'

# Install optparse
RUN apt-get install -y r-cran-optparse

# Install qualpalr from GitHub instead of CRAN
RUN R -e 'options(warn = 2); devtools::install_github("jolars/qualpalr@v.0.4.3")'

# Set R library path to avoid binding from host system
ENV R_LIBS_USER /usr/local/lib/R/site-library
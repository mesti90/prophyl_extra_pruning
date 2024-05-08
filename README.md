
<!-- README.md is generated from README.Rmd. Please edit that file -->

# prophyl

This repository is part of a set of repositories we have developed for
analysing *Acinetobacter baumannii* genomes. The master repository of
the project is at
<https://github.com/stitam/Koncz-et-al-Genomic-surveillance>.

This repository contains code for analysis steps that are not *A.
baumannii* specific, e.g. 1. building time-calibrated phylogenetic trees
(build_tree.nf) 2. analysing global and regional transmission dynamics
(get_risks.nf), 3. inferring transmission events using ancestral
geographic reconstruction (get_state_changes.nf).

The philosophy of the master repository and this repository are slightly
different: while the aim for the master repository is *only*
reproducibility (for the purpose of the article), the aim for this
repository is reusability.

At its current state the pipelines in this repository should work with
any set of prokaryotic genomes. Ongoing work includes tidying the
repository to ensure it can be run easily on other platforms.

## Repository Structure

This repository contains:

- A pipeline launcher for building time-calibrated phylogenetic trees
  (build_tree.nf)
- A pipeline launcher for analysing global and regional transmission
  dynamics (get_risks.nf)
- A pipeline launcher for inferring transmission events using ancestral
  geographic reconstruction (get_state_changes.nf)

The pipelines use both published software tools and also code developed
for this project.

- Scripts used or intended to be used within any of the automated
  pipelines are stored in `bin`.
- Other scripts are stored in `scripts`.
- R functions that are used in R scripts are stored in `R`, their man
  pages are stored in `man`.

## Software Dependencies

The automated pipelines run under Linux and require Java 11 or later as
well as Singularity (<https://sylabs.io/singularity/>). Processes in the
pipelines are fully containerised. For container versions, see the
section “Container versions” in the beginning of each pipeline launcher.

R scripts and R functions are OS independent and may run on any
platform. Furthermore, R functions are organised in an R package-like
structure. Currently you can load these functions into an R session by
running:

``` r
library(devtools)
load_all("<path to downloaded repo>")
```

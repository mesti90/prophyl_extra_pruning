#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// This pipeline takes a phylogenetic tree, attempts to date and plot it. This
// pipeline is temporary and will be integrated into another pipeline.

// Inputs - define these in a yaml file in each run
// assembly_summary = "aci_all.rds"
// tree = "tree.nwk"
// snps = "snpsites.fasta"
// poppunk = "ppdb_clusters.csv"

// Containers
r_container = "stitam/r-aci:0.1"

process date_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/post/date_tree"

    output:
    path "treedater_tree_with_time.nwk"
    path("dated_tree.rds")
    path("treedater_log.txt")
    path("treedater_root_to_tip.pdf")
    path("treedater_root_to_tip.png")
    path("treedater_tree_with_time.nwk") 

    script:
    """
    Rscript $projectDir/bin/date_tree.R \
    $launchDir/$params.tree \
    $launchDir/$params.snps \
    $launchDir/$params.assembly_summary \
    ${task.cpus}
    """
}

process plot_tree_fan {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/post/plot_tree_fan"

    input:
    path tree

    output:
    path "tree_all_tips.pdf"
    path "tree_dropped_tips.pdf"

    script:
    """
    Rscript $projectDir/bin/plot_tree_fan.R \
    $projectDir \
    $tree \
    $params.assembly_summary \
    $params.poppunk
    """
}

workflow {
    // date tree
    date_tree()
    // plot dated tree
    date_tree.out[0] | plot_tree_fan
}
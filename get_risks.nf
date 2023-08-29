#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// CONTAINERS
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/prophyl:0.10"

// PARAMETERS

// Input files
params.assemblies = "${launchDir}/assemblies.tsv"
params.tree = "${launchDir}/results/shrink_tree/treeshrink.tre"
params.snps = "${launchDir}/results/build_tree/chromosomes.nodup.filtered_polymorphic_sites.fasta"
params.dated_tree = "${launchDir}/results/choose_dated_tree/final_dated_tree.rds"

// Number and size of subsampled trees
params.subsample_count = 25
params.subsample_tipcount = 500

// Number of trees to simulate from each subsampled tree
// Simulated trees will have the same topology but different branch lengths
// Branch lengths will be simulated based on tree dating results
params.simtrees = 1

// Number of bootstrap events to perform on each tree
// Used for calculating relative risks
params.nboot_on_simtree = 1

// Variable in the assembly table to focus risk analysis on
params.focus_by = "continent"

// Value of the variable chosen above to focus on
params.focus_on = "europe"

// Processes 

process build_subset_tree {
    container "$fasttree_container"
    storeDir "$launchDir/results/build_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps)

    output:
    tuple val(subset_id), path(subset_snps), path("${subset_id}.nwk")         

    script:
    """
    FastTree -nt $subset_snps > "${subset_id}.nwk"
    """
}

process calculate_distances {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/calculate_distances"

    input:
    path simtree_paths

    output:
    tuple path("same_country.rds"),
          path("neighbors.rds"),
          path("same_continent.rds"),
          path("geodist.rds"),
          path("colldist.rds"),
          path("phylodist_list.rds")

    script:
    """
    Rscript $projectDir/bin/calculate_distances.R \
    $projectDir \
    $params.assemblies \
    $simtree_paths \
    $params.focus_by \
    $params.focus_on
    """
}

process calculate_relative_risks {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/calculate_relative_risks"

    input:
    tuple path(same_country), path(neighbors), path(same_continent), path(geodist), path(colldist), path(phylodist_list)

    output:
    path "relative_risks.rds"
    path "relative_risks.pdf"
    path "relative_risks.png"

    script:
    """
    Rscript $projectDir/bin/calculate_relative_risks.R \
    $params.assemblies \
    $colldist \
    $same_country \
    $same_continent \
    $geodist \
    $phylodist_list \
    $params.nboot_on_simtree \
    $projectDir
    """
}

process choose_dated_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/choose_dated_subset_tree/${subset_id}"

    input:
    tuple val(subset_id), path(dated_trees)

    output:
    tuple val(subset_id), path("final_dated_tree.rds"), emit: dated_tree
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/choose_dated_tree.R --trees $dated_trees
    """
}

process date_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/date_subset_tree/${subset_id}"

    input:
    tuple val(subset_id), path(subset_snps), path(subset_trees)

    output:
    tuple val(subset_id), path("dated_trees.rds"), emit: dated_trees
    path "rtt_plots/*.pdf"
    path "dated_trees/*.tre"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/date_tree.R \
    --project_dir $projectDir \
    --trees $subset_trees \
    --snps $subset_snps \
    --assemblies $params.assemblies \
    --threads ${task.cpus} \
    --branch_dimension snp_per_genome \
    --reroot false
    """
}

process filter_snps {
    //TODO create container from scratch
    container "staphb/snp-sites:2.5.1"
    storeDir "$launchDir/results/filter_snps"

    input:
    tuple val(subsample_id), path(alignment)

    output:
    tuple val(subsample_id), path("${subsample_id}_filtered.fasta")

    script:
    """
    snp-sites -o "${subsample_id}_filtered.fasta" $alignment
    """
}

process root_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/root_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps), path(subset_tree)

    output:
    tuple val(subset_id), path(subset_snps), path("rooted_trees_${subset_id}.rds"), emit: rooted_trees
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/root_subset_tree.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --dated_tree $params.dated_tree \
    --subset_tree $subset_tree \
    --threads ${task.cpus}
    """
}

process simulate_subset_trees {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/simulate_subset_trees"

    input:
    tuple val(subset_id), path(subset_tree_rds)

    output:
    path "${subset_id}.txt"
    path "${subset_id}.rds"

    script:
    """
    Rscript $projectDir/bin/simulate_subset_trees.R \
    $subset_id \
    $subset_tree_rds \
    $params.simtrees \
    ${task.cpus} \
    $launchDir
    """
}

process subsample_input {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/subsample_input"

    input:
    path assemblies

    output:
    path "subsample_*.tsv"

    script:
    """
    Rscript $projectDir/bin/subsample_input.R \
    $assemblies \
    $params.tree \
    $params.subsample_count \
    $params.subsample_tipcount
    """
}

process subset_snps {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/subset_snps"

    input:
    tuple val(subsample_id), path(subsample)

    output:
    tuple val(subsample_id), path("${subsample_id}.fasta")

    script:
    """
    Rscript $projectDir/bin/subset_snps.R \
    $params.snps \
    $subsample
    """
}

process validate_input {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/validate_input"

    output:
    path "assemblies.tsv"

    script:
    """
    Rscript $projectDir/bin/validate_input.R $params.assemblies
    """
}

// Workflow

workflow {
    validate_input()
    // prepare random subsamples from assemblies
    validate_input.out | subsample_input
    // create a channel from random subsamples
    subsample_ch = subsample_input.out.flatten() | map { [it.getBaseName(), it] }
    // build subset trees, root them
    subsample_ch | subset_snps | filter_snps | build_subset_tree | root_subset_tree 
    // date rooted subset trees
    root_subset_tree.out.rooted_trees | date_subset_tree
    // choose a single dated tree for each subset
    date_subset_tree.out.dated_trees | choose_dated_subset_tree
    // simulate trees from each subset tree
    choose_dated_subset_tree.out.dated_tree |simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simtree_paths = simulate_subset_trees.out[0].collectFile(
        name: "simtree_paths.txt",
        storeDir: "$launchDir/results/"
    )
    simtree_paths | calculate_distances | calculate_relative_risks
}

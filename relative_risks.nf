#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Containers
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/r-prophyl:0.6"

// Parameters
params.tree = "${launchDir}/results/shrink_snp_cols/treeshrink.tre"
params.snps = "${launchDir}/results/shrink_snp_cols/shrinked_snp_cols.fasta"

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

process collapse_outbreaks {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/collapse_outbreaks"

    input:
    path assemblies

    output:
    path "assemblies_collapsed_outbreaks.rds"

    script:
    """
    Rscript $projectDir/bin/collapse_outbreaks.R \
    $assemblies \
    geo_date
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

process date_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/date_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps), path(subset_tree)

    output:
    tuple val(subset_id),
          path("${subset_id}/dated_tree.rds"), \
          path("${subset_id}/treedater_log.txt"), \
          path("${subset_id}/treedater_root_to_tip.pdf"), \
          path("${subset_id}/treedater_root_to_tip.png"), \
          path("${subset_id}/treedater_tree_with_time.nwk") 

    script:
    """
    Rscript $projectDir/bin/date_subset_tree.R $subset_id $subset_tree $subset_snps $params.assemblies ${task.cpus}
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

process simulate_subset_trees {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/simulate_subset_trees"

    input:
    tuple val(subset_id), path(subset_tree_rds), path(B), path(C), path(D), path(E)

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
    // collapse outbreaks prepare random subsamples from assemblies
    validate_input.out | collapse_outbreaks | subsample_input
    // create a channel from random subsamples
    subsample_ch = subsample_input.out.flatten() | map { [it.getBaseName(), it] }
    // build subset trees, date subset trees, simulate new trees using the dated trees
    subsample_ch | subset_snps | filter_snps | build_subset_tree | date_subset_tree | simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simtree_paths = simulate_subset_trees.out[0].collectFile(
        name: "simtree_paths.txt",
        storeDir: "$launchDir/results/"
    )
    simtree_paths | calculate_distances | calculate_relative_risks
}

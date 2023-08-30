#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// CONTAINERS
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/prophyl:0.10"

// PARAMETERS

// Input files
params.assemblies = "${launchDir}/assemblies.tsv"
params.tree = "${launchDir}/results/add_duplicates/dated_tree.rds"
params.snps = "${launchDir}/results/build_tree/chromosomes.nodup.filtered_polymorphic_sites.fasta"
params.duplicates = "${launchDir}/results/remove_duplicates/duplicates.txt"

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
// Can be either "none" in which case there will be no focus
// or a variable of the assembly table e.g. "continent"
params.focus_by = "continent"

// Value of the variable chosen above to focus on
// Can be either "none" in which case there will be no focus
// or a value of the variable chosen above e.g. "europe"
params.focus_on = "europe"

// Processes 

process add_subset_duplicates {
    container "$r_container"
    storeDir "$launchDir/results/add_subset_duplicates/${subset_id}"

    input:
    tuple val(subset_id), path(dated_tree), path(duplicates)

    output:
    tuple val(subset_id), path("dated_tree.rds")

    script:
    """
    Rscript $projectDir/bin/add_duplicates.R \
    --project_dir $projectDir \
    --tree $dated_tree \
    --duplicates $duplicates
    """
}

process build_subset_tree {
    container "$fasttree_container"
    storeDir "$launchDir/results/build_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps), path(duplicates)

    output:
    tuple val(subset_id), path(subset_snps), path("${subset_id}.nwk"), path(duplicates)      

    script:
    """
    FastTree -nt $subset_snps > "${subset_id}.nwk"
    """
}

process choose_dated_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/choose_dated_subset_tree/${subset_id}"

    input:
    tuple val(subset_id), path(dated_trees), path(duplicates)

    output:
    tuple val(subset_id), path("final_dated_tree.rds"), path(duplicates), emit: dated_tree
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
    tuple val(subset_id), path(subset_snps), path(subset_trees), path(duplicates)

    output:
    tuple val(subset_id), path("dated_trees.rds"), path(duplicates), emit: dated_trees
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
    --branch_dimension snp_per_site \
    --reroot false
    """
}

process filter_snps {
    //TODO create container from scratch
    container "staphb/snp-sites:2.5.1"
    storeDir "$launchDir/results/filter_snps"

    input:
    tuple val(subsample_id), path(alignment), path(duplicates)

    output:
    tuple val(subsample_id), path("${subsample_id}_filtered.fasta"), path(duplicates)

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
    tuple val(subset_id), path(subset_snps), path(subset_tree), path(duplicates)

    output:
    tuple val(subset_id), path(subset_snps), path("rooted_trees_${subset_id}.rds"), path(duplicates), emit: rooted_trees
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/root_subset_tree.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --dated_tree $params.tree \
    --subset_tree $subset_tree \
    --threads ${task.cpus}
    """
}

process rr_calc_counts {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/rr_calc_counts"

    input:
    tuple \
      path(same_city),
      path(same_country),
      path(neighbors),
      path(same_continent),
      path(geodist),
      path(colldist),
      path(phylodist_list)

    output:
    path "countlist.rds"

    script:
    """
    Rscript $projectDir/bin/rr_calc_counts.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --colldist $colldist \
    --same_city $same_city \
    --same_country $same_country \
    --neighbors $neighbors \
    --same_continent $same_continent \
    --geodist $geodist \
    --phylodist $phylodist_list \
    --nboot $params.nboot_on_simtree
    """
}

process rr_calc_dist {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/rr_calc_dist"

    input:
    path simtree_paths

    output:
    tuple \
      path("same_city.rds"),
      path("same_country.rds"),
      path("neighbors.rds"),
      path("same_continent.rds"),
      path("geodist.rds"),
      path("colldist.rds"),
      path("phylodist_list.rds")

    script:
    """
    Rscript $projectDir/bin/rr_calc_dist.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --simtrees $simtree_paths \
    --focus_by $params.focus_by \
    --focus_on $params.focus_on
    """
}

process rr_plot_risks {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/rr_plot_risks"

    input:
    path countlist

    output:
    path "relative_risks_type1.rds"
    path "relative_risks_type1.pdf"
    path "relative_risks_type1.png"
    path "relative_risks_type2.rds"
    path "relative_risks_type2.pdf"
    path "relative_risks_type2.png"
    path "relative_risks_type3.rds"
    path "relative_risks_type3.pdf"
    path "relative_risks_type3.png"

    script:
    """
    Rscript $projectDir/bin/rr_plot_risks.R --countlist $countlist
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
    path "subsample_*.tsv", emit: subsample
    path "subsample_*.rds", emit: duplicates

    script:
    """
    Rscript $projectDir/bin/subsample_input.R \
    --project_dir $projectDir \
    --assemblies $assemblies \
    --tree $params.tree \
    --duplicates $params.duplicates \
    --subsample_count $params.subsample_count \
    --subsample_tipcount $params.subsample_tipcount
    """
}

process subset_snps {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/subset_snps"

    input:
    tuple val(subsample_id), path(subsample), path(duplicates)

    output:
    tuple val(subsample_id), path("${subsample_id}.fasta"), path(duplicates)

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
    subsample_ch = subsample_input.out.subsample.flatten() | map { [it.getBaseName(), it] }
    duplicate_ch = subsample_input.out.duplicates.flatten() | map { [it.getBaseName(), it] }
    subsample_tuple_ch = subsample_ch.join(duplicate_ch)
    // build subset trees, root them
    subsample_tuple_ch | subset_snps | filter_snps | build_subset_tree | root_subset_tree 
    // date rooted subset trees
    root_subset_tree.out.rooted_trees | date_subset_tree
    // choose a single dated tree for each subset
    date_subset_tree.out.dated_trees | choose_dated_subset_tree
    // simulate trees from each subset tree
    choose_dated_subset_tree.out.dated_tree |add_subset_duplicates | simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simtree_paths = simulate_subset_trees.out[0].collectFile(
        name: "simtree_paths.txt",
        storeDir: "$launchDir/results/"
    )
    simtree_paths | rr_calc_dist | rr_calc_counts | rr_plot_risks
}

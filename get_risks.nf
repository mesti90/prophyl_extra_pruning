#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// CONTAINERS
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/prophyl:0.11"

// PARAMETERS

// Output parameters

params.resdir = "results"

// Input files
params.ASSEMBLIES = "${launchDir}/assemblies.tsv"
params.assemblies = "${launchDir}/assemblies_for_country_rr.tsv"
params.tree = "${launchDir}/$params.resdir/add_duplicates/dated_tree.rds"
params.snps = "${launchDir}/$params.resdir/build_tree/chromosomes.nodup.filtered_polymorphic_sites.fasta"
params.duplicates = "${launchDir}/$params.resdir/remove_duplicates/duplicates.txt"

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
params.focus_by = "none"

// Value of the variable chosen above to focus on
// Can be either "none" in which case there will be no focus
// or a value of the variable chosen above e.g. "europe"
params.focus_on = "none"

// MRCA categories for relative risk analysis
params.mrca_categories = "0,6,12,40"

// Maximum distance in collection dates between two isolates
// to be considered as a pair for risk analysis
params.colldist_max = 2

// Processes 

process add_subset_duplicates {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/add_subset_duplicates/${subset_id}"

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
    storeDir "$launchDir/$params.resdir/build_subset_tree"

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
    storeDir "$launchDir/$params.resdir/choose_dated_subset_tree/${subset_id}"

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
    storeDir "$launchDir/$params.resdir/date_subset_tree/${subset_id}"

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
    --assemblies $params.ASSEMBLIES \
    --threads ${task.cpus} \
    --branch_dimension snp_per_site \
    --reroot false
    """
}

process filter_snps {
    //TODO create container from scratch
    container "staphb/snp-sites:2.5.1"
    storeDir "$launchDir/$params.resdir/filter_snps"

    input:
    tuple val(subsample_id), path(alignment), path(duplicates)

    output:
    tuple val(subsample_id), path("${subsample_id}_filtered.fasta"), path(duplicates)

    script:
    """
    snp-sites -o "${subsample_id}_filtered.fasta" $alignment
    """
}

process qc_dated_subset_tree {
    container "$r_container"                                                                                                                                                                                                                                                                                                            
    storeDir "$launchDir/$params.resdir/qc_dated_subset_tree/${subset_id}"

    input:
    tuple val(subset_id), path(dated_tree), path(duplicates)

    output:
    path "${subset_id}_qc.tsv"

    script:
    """
    Rscript $projectDir/bin/qc_dated_subset_tree.R \
    --project_dir $projectDir \
    --subsample_id $subset_id \
    --tree $dated_tree
    """
}

process root_subset_tree {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/root_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps), path(subset_tree), path(duplicates)

    output:
    tuple val(subset_id), path(subset_snps), path("rooted_trees_${subset_id}.rds"), path(duplicates), emit: rooted_trees
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/root_subset_tree.R \
    --project_dir $projectDir \
    --assemblies $params.ASSEMBLIES \
    --dated_tree $params.tree \
    --subset_tree $subset_tree \
    --threads ${task.cpus}
    """
}

process rr_calc_counts {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/rr_calc_counts"

    input:
    tuple \
      path(same_city),
      path(same_country),
      path(neighbors),
      path(same_continent),
      path(geodist),
      path(colldist),
      path(phylodist),
      path(phylodist_list)

    output:
    tuple \
      path("countlist_all.rds"),
      path("countlist.rds") 

    script:
    """
    Rscript $projectDir/bin/rr_calc_counts.R \
    --project_dir $projectDir \
    --assemblies $params.ASSEMBLIES \
    --assemblies_collapsed $params.assemblies \
    --colldist $colldist \
    --same_city $same_city \
    --same_country $same_country \
    --neighbors $neighbors \
    --same_continent $same_continent \
    --geodist $geodist \
    --phylodist_all $phylodist \
    --phylodist $phylodist_list \
    --nboot $params.nboot_on_simtree \
    --mrca_categories $params.mrca_categories \
    --colldist_max $params.colldist_max
    """
}

process rr_calc_dist {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/rr_calc_dist"

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
      path("phylodist.rds"),
      path("phylodist_list.rds")

    script:
    """
    Rscript $projectDir/bin/rr_calc_dist.R \
    --project_dir $projectDir \
    --assemblies $params.ASSEMBLIES \
    --tree $params.tree \
    --simtrees $simtree_paths \
    --focus_by $params.focus_by \
    --focus_on $params.focus_on
    """
}

process rr_plot_risks {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/rr_plot_risks"

    input:
    tuple path(countlist_all), path(countlist)

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
    path "p_values.tsv"

    script:
    """
    Rscript $projectDir/bin/rr_plot_risks.R --countlist $countlist --countlist_all $countlist_all
    """
}

process simulate_subset_trees {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/simulate_subset_trees"

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
    $launchDir \
    $launchDir/$params.resdir
    """
}

process subsample_input {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/subsample_input"

    input:
    path assemblies

    output:
    path "subsample_*.tsv", emit: subsample
    path "subsample_*.rds", emit: duplicates

    script:
    """
    Rscript $projectDir/bin/subsample_input.R \
    --project_dir $projectDir \
    --ASSEMBLIES $params.ASSEMBLIES \
    --assemblies $assemblies \
    --tree $params.tree \
    --duplicates $params.duplicates \
    --subsample_count $params.subsample_count \
    --subsample_tipcount $params.subsample_tipcount \
    --focus_by $params.focus_by \
    --focus_on $params.focus_on
    """
}

process subset_snps {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/subset_snps"

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
    storeDir "$launchDir/$params.resdir/validate_input"

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
    // get QC metrics for dated trees
    choose_dated_subset_tree.out.dated_tree | qc_dated_subset_tree
    qc_dated_subset_tree.out.collectFile(
        name: "qc_dated_trees.tsv",
        storeDir: "$launchDir/$params.resdir/",
        keepHeader: true
    )
    // simulate trees from each subset tree
    choose_dated_subset_tree.out.dated_tree |add_subset_duplicates | simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simtree_paths = simulate_subset_trees.out[0].collectFile(
        name: "simtree_paths.txt",
        storeDir: "$launchDir/$params.resdir/"
    )
    simtree_paths | rr_calc_dist | rr_calc_counts | rr_plot_risks
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// assemblies_ch = Channel.fromPath("$launchDir/assemblies.tsv")

ans_targets = Channel.of("country", "continent", "mlst", "k_serotype")

// tree_ch = Channel.fromPath("$launchDir/treedater_tree_with_time.nwk", checkIfExists: true)
// pastml_ch = Channel.fromPath("$launchDir/output/pastml/combined_ancestral_states.tab", checkIfExists: true)
// meta_rda_ch = Channel.fromPath("$launchDir/treemeta.rda", checkIfExists: true)
// meta_tsv_ch = Channel.fromPath("$launchDir/treemeta.tsv", checkIfExists: true)

process bootstrap_tree {
    container "staphb/iqtree"
    storeDir "$launchDir/results/bootstrap_tree"

    input:
    tuple path(A), path(B), path(snps), path(D), path(tree), path(F), path(G), path(H), path(I) 

    output:
    path "chromosomes.filtered_polymorphic_sites.fasta.treefile"

    script:
    """
    iqtree \
    -t $tree \
    -s $snps \
    -nt ${task.cpus} \
    -mem "${task.memory.toGiga()}G" \
    -bb $params.bootstrap_replicates \
    -wbtl
    """
}

process build_tree {
    container "mesti90/gubbins:latest"
    storeDir "$launchDir/results/build_tree"

    input:
    path chromosomes

    output:
    tuple path("chromosomes.log"), \
          path("chromosomes.branch_base_reconstruction.embl"), \
          path("chromosomes.filtered_polymorphic_sites.fasta"), \
          path("chromosomes.filtered_polymorphic_sites.phylip"), \
          path("chromosomes.node_labelled.final_tree.tre"), \
          path("chromosomes.per_branch_statistics.csv"), \
          path("chromosomes.recombination_predictions.embl"), \
          path("chromosomes.recombination_predictions.gff"), \
          path("chromosomes.summary_of_snp_distribution.vcf")          

    script:
    """
    run_gubbins.py \
    --model-fitter raxml \
    --tree-builder fasttree \
    --threads ${task.cpus} \
    --iterations $params.gubbins_iterations\
    $chromosomes
    """
}

process build_subset_tree {
    container "staphb/fasttree:latest"
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
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/calculate_distances"

    input:
    path simtree_paths

    output:
    tuple path("same_country.rds"), path("neighbors.rds"), path("same_continent.rds"), path("geodist.rds"), path("colldist.rds"), path("phylodist_list.rds")

    script:
    """
    Rscript $projectDir/bin/calculate_distances.R $params.Rdir $projectDir $params.assemblies $simtree_paths
    """
}

process calculate_relative_risks {
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/calculate_relative_risks"

    input:
    tuple path(same_country), path(neighbors), path(same_continent), path(geodist), path(colldist), path(phylodist_list)

    output:
    path "risk_ratios.rds"

    script:
    """
    Rscript $projectDir/bin/calculate_relative_risks.R $params.assemblies $colldist $same_country $same_continent $geodist $phylodist_list $params.nboot_on_simtree
    """
}

process create_genome_list {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/create_genome_list"

    input:
    path assemblies

    output:
    file "log.txt"
    file "paired_reads.csv"
    file "single_reads.csv"
    file "contigs.csv"

    script:
    """
    Rscript $projectDir/bin/prep_snippy_input.R $assemblies "$launchDir/genomes"
    """
}

process date_tree {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/date_tree"

    input:
    tuple path(A), path(B), path(snps), path(D), path(tree), path(F), path(G), path(H), path(I) 

    output:
    path "treedater_tree_with_time.nwk"
    tuple path("dated_tree.rds"), \
          path("treedater_log.txt"), \
          path("treedater_root_to_tip.pdf"), \
          path("treedater_root_to_tip.png"), \
          path("treedater_tree_with_time.nwk") 

    script:
    """
    Rscript $projectDir/bin/date_tree.R $tree $snps $params.assemblies ${task.cpus}
    """
}

process date_subset_tree {
    container "stitam/r-packages:1.8"
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

process keep_chromosome {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/keep_chromosome"

    input:
    path assembly_dir

    output:
    path "${assembly_dir}.fasta"

    script:
    """
    Rscript $projectDir/bin/keep_chromosome.R $assembly_dir
    """
}

process plot_tree {
    container "stitam/rplots:1.0"
    storeDir "$launchDir/results/plot_tree"

    input:
    path tree_tbl

    output:
    path "dated_tree.pdf"
    path "dated_tree.png"
    
    script:
    """
    Rscript $projectDir/bin/plot_tree.R $tree_tbl
    """
}


process predict_ancestral_states {
    container "evolbioinfo/pastml"
    storeDir "$launchDir/results/predict_ancestral_states"
    
    input:
    tuple path(dated_tree), val(target)

    output:
    tuple val(target), path(target)

    script:
    """
    pastml \
    -t $dated_tree \
    -d $params.assemblies \
    -c $target \
    --threads ${task.cpus} \
    --work_dir $target \
    --offline
    """
}

process prep_tree_tbl {
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/prep_tree_tbl"

    input:
    path tree
    path ancestral_states

    output:
    path "tree_tbl.rds"
    
    script:
    """
    Rscript $projectDir/bin/prep_tree_tbl.R $tree $params.assemblies $ancestral_states $params.Rdir
    """
}

process simulate_trees {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/simulate_trees"

    input:
    tuple path(tree_rds), path(B), path(C), path(D), path(E)

    output:
    path "simtrees.rds"

    script:
    """
    Rscript $projectDir/bin/simulate_trees.R $tree_rds $params.simtrees ${task.cpus}
    """
}

process simulate_subset_trees {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/simulate_subset_trees"

    input:
    tuple val(subset_id), path(subset_tree_rds), path(B), path(C), path(D), path(E)

    output:
    path "${subset_id}.txt"
    path "${subset_id}.rds"

    script:
    """
    Rscript $projectDir/bin/simulate_subset_trees.R $subset_id $subset_tree_rds $params.simtrees ${task.cpus} $launchDir
    """
}

process snippy_contig {
    container "staphb/snippy"
    storeDir "$launchDir/results/snippy_contig"

    input:
    tuple val(assembly_id), val(contigs)

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $params.reffile \
    --ctgs $contigs \
    --force
    """
}

process snippy_paired {
    container "staphb/snippy"
    storeDir "$launchDir/results/snippy_paired"

    input:
    tuple val(assembly_id), val(R1), val(R2)

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $params.reffile \
    --R1 $R1 \
    --R2 $R2 \
    --force
    """
}

process snippy_single {
    container "staphb/snippy"
    storeDir "$launchDir/results/snippy_single"

    input:
    tuple val(assembly_id), val(reads)

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $params.reffile \
    --se $reads \
    --force
    """
}

process subsample_input {
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/subsample_input"

    input:
    path assemblies

    output:
    path "subsample_*.tsv"

    script:
    """
    Rscript $projectDir/bin/subsample_input.R $assemblies $params.subsample_count $params.subsample_tipcount
    """
}

process subset_snps {
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/subset_snps"

    input:
    tuple path(A), path(B), path(snps), path(D), path(E), path(F), path(G), path(H), path(I), val(subsample_id), path(subsample)

    output:
    tuple val(subsample_id), path("${subsample_id}.fasta")

    script:
    """
    Rscript $projectDir/bin/subset_snps.R $snps $subsample $params.Rdir
    """
}

process tidy_ancestral_states {
    container "stitam/r-packages:1.8"
    storeDir "$launchDir/results/tidy_ancestral_states"

    input:
    tuple val(target), path(combined)

    output:
    path "${target}.tsv"

    script:
    """
    Rscript $projectDir/bin/tidy_ancestral_states.R $combined $target
    """
}

process tidy_bootstrap_tree {
    container "stitam/r-bio:1.0"
    storeDir "$launchDir/results/tidy_bootstrap_tree"

    input: 
    path bstree

    output:
    path "tree_tbl.rds"

    script:
    """
    Rscript $projectDir/bin/tidy_bootstrap_tree.R $bstree $params.Rdir
    """
} 

process validate_input {
    container "stitam/r-bio:1.1" 
    storeDir "$launchDir/results/validate_input"

    output:
    path "assemblies.tsv"

    script:
    """
    Rscript $projectDir/bin/validate_input.R $params.assemblies
    """
}

workflow {
    // Prepare pseudo-whole genomes
    validate_input() | create_genome_list
    create_genome_list.out[1].splitCsv(header: true) | snippy_paired
    create_genome_list.out[2].splitCsv(header: true) | snippy_single
    create_genome_list.out[3].splitCsv(header: true) | snippy_contig
    snippy_paired.out.concat(snippy_single.out, snippy_contig.out) | keep_chromosome
    // Mask recombination, build tree, bootstrap
    keep_chromosome.out.collectFile(name: "chromosomes.fasta") | build_tree | bootstrap_tree | tidy_bootstrap_tree
    // Date tree with treedater
    build_tree.out | date_tree
    // Simulate new trees using the dated tree, calculate geo distance and phylo distance, calculate relative_risks
    // date_tree.out[1] | simulate_trees | calculate_distances | calculate_relative_risks
    // predict ancestral states for each variable defined in the ans_targets channel
    date_tree.out[0].combine(ans_targets) | predict_ancestral_states | tidy_ancestral_states
    // prepare tree_tbl which contains predicted ancestral states for internal nodes
    prep_tree_tbl(date_tree.out[0], tidy_ancestral_states.out.collectFile(name: "all_ancestral_states.tsv", newLine: false, keepHeader: true))
    // prepare random subsamples from assemblies and create a channel
    validate_input.out | subsample_input
    subsample_ch = subsample_input.out.flatten() | map { [it.getBaseName(), it] }
    // build subset trees, date subset trees, simulate new trees using the dated trees
    build_tree.out.combine(subsample_ch) | subset_snps | build_subset_tree | date_subset_tree | simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simulate_subset_trees.out[0].collectFile(name: "simtree_paths.txt") | calculate_distances | calculate_relative_risks
}
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// assemblies_ch = Channel.fromPath("$launchDir/assemblies.tsv")

ans_targets = Channel.of("country", "continent", "mlst", "k_serotype")

// tree_ch = Channel.fromPath("$launchDir/treedater_tree_with_time.nwk", checkIfExists: true)
// pastml_ch = Channel.fromPath("$launchDir/output/pastml/combined_ancestral_states.tab", checkIfExists: true)
// meta_rda_ch = Channel.fromPath("$launchDir/treemeta.rda", checkIfExists: true)
// meta_tsv_ch = Channel.fromPath("$launchDir/treemeta.tsv", checkIfExists: true)

// Containers
gubbins_container = "mesti90/gubbins:latest"
hgttree_container = "mesti90/hgttree:2.11"
iqtree_container = "staphb/iqtree"
fasttree_container = "staphb/fasttree:latest"
pastml = "evolbioinfo/pastml"
r_container = "stitam/r-prophyl:0.7"
root_digger_container = "stitam/root_digger:1.7.0"
snippy_container = "staphb/snippy"

process bootstrap_tree {
    container "$iqtree_container"
    storeDir "$launchDir/results/bootstrap_tree"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D)  

    output:
    path "chromosomes.filtered_polymorphic_sites.fasta.treefile"

    script:
    """
    iqtree \
    -t $shrinked_tree \
    -s $shrinked_snps \
    -nt ${task.cpus} \
    -mem "${task.memory.toGiga()}G" \
    -bb $params.bootstrap_replicates \
    -wbtl
    """
}

process build_subset_tree {
    container "$fasttree_container"
    storeDir "$launchDir/results/build_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps)

    output:
    tuple val(subset_id), path(subset_snps), path("${subset_id}.nwk"), emit: subset_tree         

    script:
    """
    FastTree -nt $subset_snps > "${subset_id}.nwk"
    """
}

process build_tree {
    container "$gubbins_container"
    storeDir "$launchDir/results/build_tree"

    input:
    path chromosomes

    output:
    tuple path("chromosomes.nodup.log"), \
          path("chromosomes.nodup.branch_base_reconstruction.embl"), \
          path("chromosomes.nodup.filtered_polymorphic_sites.fasta"), \
          path("chromosomes.nodup.filtered_polymorphic_sites.phylip"), \
          path("chromosomes.nodup.node_labelled.final_tree.tre"), \
          path("chromosomes.nodup.per_branch_statistics.csv"), \
          path("chromosomes.nodup.recombination_predictions.embl"), \
          path("chromosomes.nodup.recombination_predictions.gff"), \
          path("chromosomes.nodup.summary_of_snp_distribution.vcf")          

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
    $projectDir \
    $neighbors
    """
}

process choose_dated_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/choose_dated_tree"

    input:
    path dated_trees

    output:
    path "final_dated_tree.rds", emit: dated_big_tree
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/choose_dated_tree.R $dated_trees
    """
}

process choose_reference_genome {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/choose_reference_genome"

    input:
    path assemblies

    output:
    path "*.*"

    script:
    """
    Rscript $projectDir/bin/choose_reference_genome.R \
    --assemblies $assemblies \
    --genome_dir $params.genome_dir
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

process create_genome_list {
    container "$r_container"
    containerOptions "--no-home"
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

process date_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/date_subset_tree"

    input:
    tuple val(subset_id), path(subset_snps), path(rooted_subset_tree)

    output:
    tuple val(subset_id),
          path("${subset_id}/dated_tree.rds"), \
          path("${subset_id}/treedater_log.txt"), \
          path("${subset_id}/treedater_root_to_tip.pdf"), \
          path("${subset_id}/treedater_root_to_tip.png"), \
          path("${subset_id}/treedater_tree_with_time.nwk") 

    script:
    """
    Rscript $projectDir/bin/date_subset_tree.R \
    $subset_id \
    $rooted_subset_tree \
    $subset_snps \
    $params.assemblies \
    ${task.cpus}
    """
}

process date_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/date_tree/$root_method"

    input:
    tuple val(root_method), path(rooted_tree), path(snps)

    output:
    path "treedater_tree_with_time.nwk"
    path "dated_tree_${root_method}.rds", emit: dated_tree_rds
    path "treedater_log.txt"
    path "treedater_root_to_tip.pdf"
    path "treedater_root_to_tip.png"

    script:
    """
    Rscript $projectDir/bin/date_tree.R \
    $rooted_tree \
    $snps \
    $params.assemblies \
    ${task.cpus} \
    snp_per_genome \
    $params.reroot_tree \
    $root_method
    """
}

process date_tree_bactdating {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/date_tree_bactdating"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D)   

    output:
    path "treedater_tree_with_time.nwk"
    tuple path("dated_tree.rds"), \
          path("dated_tree.tre"), \
          path("trace.pdf"), \
          path("root_to_tip_regression.pdf"), \
          path("log.txt")

    script:
    """
    Rscript $projectDir/bin/date_tree_bactdating.R \
    $projectDir \
    $shrinked_tree \
    $shrinked_snps \
    $params.assemblies \
    snp_per_genome
    """
}

process filter_input {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/filter_input"

    input:
    path assembly_tbl

    output:
    path "filtered_assemblies.rds", emit: filtered_assemblies
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/filter_input.R --assembly_file $assembly_tbl
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

process keep_chromosome {
    container "$r_container"
    containerOptions "--no-home"
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
    container "$r_container"
    containerOptions "--no-home"
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
    container "$pastml_container"
    containerOptions "--no-home"
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
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/prep_tree_tbl"

    input:
    path tree
    path ancestral_states

    output:
    path "tree_tbl.rds"
    
    script:
    """
    Rscript $projectDir/bin/prep_tree_tbl.R $tree $params.assemblies $ancestral_states
    """
}

process remove_duplicates {
    container "$hgttree_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/remove_duplicates"
    
    input:
    path chromosomes
    
    output:
    path "chromosomes.nodup.fasta", emit: chromosomes_nodup
    path "duplicates.txt"
    
    script:
    """
    touch duplicates.txt
    seqkit rmdup -s $chromosomes -D duplicates.txt -o chromosomes.nodup.fasta -j ${task.cpus}
    """
}

process root_subset_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/root_subset_tree/${subset_id}"

    input:
    path dated_big_tree
    tuple val(subset_id), path(subset_snps), path(subset_tree)

    output:
    path "log.txt"
    tuple \
      val(subset_id), \
      path(subset_snps), \
      path("rooted_subset_tree.rds"), \
      emit: rooted_subset_tree

    script:
    """
    Rscript $projectDir/bin/root_subset_tree.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --dated_tree $dated_big_tree \
    --subset_tree $subset_tree \
    --threads ${task.cpus}
    """
}

mad_ch = Channel.of("mad")
process root_tree_mad {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/root_tree_mad"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D)
    val root_method

    output:
    tuple val(root_method), path("rooted_tree_mad.tre"), path(shrinked_snps)
    path("rooted_tree_mad.rds")
    path("distmat_t.rds")
    path("distmat_t2.rds")
    path("log.txt")

    script:
    """
    Rscript $projectDir/bin/root_tree_mad.R \
    $projectDir \
    $shrinked_tree \
    $task.cpus
    """
}

mp_ch = Channel.of("midpoint")
process root_tree_mp {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/root_tree_mp"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D)
    val root_method

    output:
    tuple val(root_method), path("rooted_tree_mp.tre"), path(shrinked_snps)
    path("rooted_tree_mp.rds")
    path("log.txt")

    script:
    """
    Rscript $projectDir/bin/root_tree_mp.R $projectDir $shrinked_tree
    """
}

process root_tree_rd {
    container "$root_digger_container"
    storeDir "$launchDir/results/root_tree"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D) 

    output:
    tuple path(shrinked_snps), \
          path("rooted_tree.rooted.tree"), \
          path("rooted_tree.ckp")

    script:
    """
    rd \
    --msa $shrinked_snps \
    --tree $shrinked_tree \
    --prefix rooted_tree \
    --threads ${task.cpus} \
    --seed 0
    """
}

process root_tree_rtt {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/root_tree_rtt"

    input:
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D) 

    output:
    path("rooted_trees/*.tre")
    path("rooted_trees.rds")
    path("rtt_metrics.rds")
    path("rtt_plots.pdf")
    path("log.txt")

    script:
    """
    Rscript $projectDir/bin/root_tree_rtt.R \
    $projectDir \
    $shrinked_tree \
    $params.assemblies \
    $task.cpus
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

process shrink_tree {
    container "$hgttree_container"
    storeDir "$launchDir/results/shrink_tree"

    input:
    tuple path(A),
          path(B),
          path(snps),
          path(D),
          path(tree),
          path(F),
          path(G),
          path(H),
          path(I) 

    output:
    tuple path(snps),
          path("treeshrink.tre"),
          path("treeshrink.txt"),
          path("treeshrink_summary.txt")    

    script:
    """
    run_treeshrink.py --tree $tree --outprefix treeshrink --force --outdir .
    """
}

process shrink_snp_rows {
    container "$r_container"
    storeDir "$launchDir/results/shrink_snp_rows"

    input:
    tuple path(snps), path(shrinked_tree), path(C), path(D)

    output:
    tuple path("shrinked_snp_rows.fasta"), path(shrinked_tree), path(C), path(D)

    script:
    """
    Rscript $projectDir/bin/shrink_snp_rows.R $shrinked_tree $snps
    """
}

process shrink_snp_cols {
    //TODO create container from scratch
    container "staphb/snp-sites:2.5.1"
    storeDir "$launchDir/results/shrink_snp_cols"

    input:
    tuple path(shrinked_snp_rows), path(shrinked_tree), path(C), path(D)  

    output:
    tuple path("shrinked_snp_cols.fasta"), path(shrinked_tree), path(C), path(D)
    path "shrinked_snp_cols.fasta"

    script:
    """
    snp-sites -o shrinked_snp_cols.fasta $shrinked_snp_rows 
    """
}

process simulate_trees {
    container "$r_container"
    containerOptions "--no-home"
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

process snippy_contig {
    container "$snippy_container"
    storeDir "$launchDir/results/snippy_contig"

    input:
    tuple val(assembly_id), val(contigs)
    each reference_genome

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $reference_genome \
    --ctgs $contigs \
    --force
    """
}

process snippy_paired {
    container "$snippy_container"
    storeDir "$launchDir/results/snippy_paired"

    input:
    tuple val(assembly_id), val(R1), val(R2)
    each reference_genome

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $reference_genome \
    --R1 $R1 \
    --R2 $R2 \
    --force
    """
}

process snippy_single {
    container "$snippy_container"
    storeDir "$launchDir/results/snippy_single"

    input:
    tuple val(assembly_id), val(reads)
    each reference_genome

    output:
    path assembly_id
    
    script:
    """
    snippy \
    --outdir $assembly_id \
    --ref $reference_genome \
    --se $reads \
    --force
    """
}

process subsample_input {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/subsample_input"

    input:
    path assemblies
    tuple path(shrinked_snps), path(shrinked_tree), path(C), path(D)

    output:
    path "subsample_*.tsv"
    

    script:
    """
    Rscript $projectDir/bin/subsample_input.R \
    $assemblies \
    $shrinked_tree \
    $params.subsample_count \
    $params.subsample_tipcount
    """
}

process subset_snps {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/subset_snps"

    input:
    tuple path(shrinked_snps), path(B), path(C), path(D), val(subsample_id), path(subsample)

    output:
    tuple val(subsample_id), path("${subsample_id}.fasta")


    script:
    """
    Rscript $projectDir/bin/subset_snps.R \
    $shrinked_snps \
    $subsample
    """
}

process tidy_ancestral_states {
    container "$r_container"
    containerOptions "--no-home"
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
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/results/tidy_bootstrap_tree"

    input: 
    path bstree

    output:
    path "tree_tbl.rds"

    script:
    """
    Rscript $projectDir/bin/tidy_bootstrap_tree.R $bstree
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

workflow {
    // Validate input and create a list of genomes to process
    validate_input() | create_genome_list
    
    // Choose a reference genome
    if (params.reference_genome != null) {
        reference_genome = "$launchDir/$params.reference_genome"
        println("Reference genome for snippy: $reference_genome")
    } else {
        reference_genome = choose_reference_genome(validate_input.out)
    }
    
    // Construct pseudo-whole genomes

    snippy_paired(create_genome_list.out[1].splitCsv(header: true), reference_genome)
    snippy_single(create_genome_list.out[2].splitCsv(header: true), reference_genome)
    snippy_contig(create_genome_list.out[3].splitCsv(header: true), reference_genome)

    // Combine snippy outputs and keep only chromosomes
    snippy_paired.out.concat(snippy_single.out, snippy_contig.out) | keep_chromosome
    chromosomes = keep_chromosome.out.collectFile(
        name: "chromosomes.fasta",
        storeDir: "$launchDir/results/"
    )

    // Remove duplicates, mask recombination, build tree, shrink, bootstrap

    chromosomes | remove_duplicates
    remove_duplicates.out.chromosomes_nodup | build_tree | shrink_tree | shrink_snp_rows | shrink_snp_cols //| bootstrap_tree | tidy_bootstrap_tree

    // Root shrinked tree using mad, midpoint, rtt
    root_tree_mad(shrink_snp_cols.out[0], mad_ch)
    root_tree_mp(shrink_snp_cols.out[0], mp_ch)
    shrink_snp_cols.out[0] | root_tree_rtt
    
    // Combine rooted trees
    rtt_trees_ch = root_tree_rtt.out[0].flatten()
    rtt_trees_ch = rtt_trees_ch | map { [it.getBaseName().split("_tree_")[1], it] }  
    rtt_trees_ch = rtt_trees_ch.combine(shrink_snp_cols.out[1])

    rooted_trees_ch = root_tree_mad.out[0].mix(
        root_tree_mp.out[0],
        rtt_trees_ch
    )

    // Date all rooted trees
    date_tree(rooted_trees_ch)

    date_tree.out.dated_tree_rds.collect() | choose_dated_tree

    // Date shrinked tree with BactDating
    // shrink_snp_cols.out[0] | date_tree_bactdating
    // Simulate new trees using the dated tree, calculate geo distance and phylo distance, calculate relative_risks
    // date_tree.out[1] | simulate_trees | calculate_distances | calculate_relative_risks
    // predict ancestral states for each variable defined in the ans_targets channel
    // date_tree.out[0].combine(ans_targets) | predict_ancestral_states | tidy_ancestral_states
    // prepare tree_tbl which contains predicted ancestral states for internal nodes
    // prep_tree_tbl(date_tree.out[0], tidy_ancestral_states.out.collectFile(name: "all_ancestral_states.tsv", newLine: false, keepHeader: true))

    // collapse outbreaks 
    validate_input.out | filter_input
    filter_input.out.filtered_assemblies | collapse_outbreaks
    // create a channel from random subsamples prepare random subsamples from assemblies
    subsample_input(collapse_outbreaks.out, shrink_snp_cols.out[0])
    subsample_ch = subsample_input.out.flatten() | map { [it.getBaseName(), it] }
    // build subset trees
    shrink_snp_cols.out[0].combine(subsample_ch) | subset_snps | filter_snps | build_subset_tree
    // root subset trees
    root_subset_tree(
        choose_dated_tree.out.dated_big_tree, 
        build_subset_tree.out.subset_tree
    )
    // date subset trees, simulate subset trees
    root_subset_tree.out.rooted_subset_tree | date_subset_tree | simulate_subset_trees
    // calculate geo distance and phylo distance, calculate relative risks
    simtree_paths = simulate_subset_trees.out[0].collectFile(
        name: "simtree_paths.txt",
        storeDir: "$launchDir/results/"
    )
    simtree_paths | calculate_distances | calculate_relative_risks
}
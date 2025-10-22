#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Containers
gubbins_container = "stitam/prophyl:0.13"
hgttree_container = "mesti90/hgttree:2.11"
iqtree_container = "staphb/iqtree"
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/prophyl:0.13"
snippy_container = "staphb/snippy"
treepruner_container = "mesti90/treepruner:1.1"

process add_duplicates {
    container "$r_container"
    storeDir "$launchDir/$params.resdir"

    input:
    path tree
    path duplicates

    output:
    path "final_tree.rds"
    path "final_tree.nwk"
    path "final_tree.tsv"
    path "rtt_plot.png"
    path "metrics.tsv"

    script:
    """
    Rscript $projectDir/bin/add_duplicates.R \
    --project_dir $projectDir \
    --launch_dir $launchDir \
    --tree $tree \
    --duplicates $duplicates
    """
}

process bootstrap_tree {
    container "$iqtree_container"
    storeDir "$launchDir/$params.resdir/bootstrap_tree"

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

process build_tree {
    container "$gubbins_container"
    storeDir "$launchDir/$params.resdir/build_tree"

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
    nohup run_gubbins.py \
    --model-fitter raxmlng \
    --tree-builder fasttree \
    --threads ${task.cpus} \
    --iterations $params.gubbins_iterations\
    $chromosomes
    """
}

process choose_dated_tree {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/choose_dated_tree"

    input:
    tuple path(dated_tree_shrink), path(dated_tree_prune)

    output:
    path "final_dated_tree.rds", emit: dated_big_tree
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/choose_dated_tree.R --trees_shrink $dated_trees --trees_prune $dated_tree_prune
    """
}

process choose_reference_genome {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/choose_reference_genome"

    input:
    path assemblies

    output:
    path "*.*"

    script:
    """
    Rscript $projectDir/bin/choose_reference_genome.R \
    --project_dir $projectDir \
    --assemblies $assemblies
    """
}

process prep_snippy_input {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/prep_snippy_input"

    input:
    path assemblies

    output:
    path "snippy_input.tsv", emit: snippy_input
    path "unzipped_assemblies"

    script:
    """
    Rscript $projectDir/bin/prep_snippy_input.R \
    --project_dir $projectDir \
    --store_dir $launchDir/$params.resdir/prep_snippy_input \
    --assemblies $assemblies
    """
}

process date_tree {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/date_tree_${label}"

    input:
    tuple val(label), path(snps), path(rooted_trees)

    output:
    path "dated_trees.${label}.rds", emit: dated_trees
    path "rtt_plots_${label}/*.pdf"
    path "dated_trees_${label}/*.tre"
    path "log.${label}.txt"

    script:
    """
    Rscript $projectDir/bin/date_tree.R \
    --project_dir $projectDir \
    --trees $rooted_trees \
    --snps $snps \
    --assemblies $params.assemblies \
    --threads ${task.cpus} \
    --branch_dimension snp_per_genome \
    --clock $params.clock \
    --reroot false
    
    mv dated_trees.rds dated_trees.${label}.rds
    mkdir -p dated_trees_${label} && mv dated_trees/* dated_trees_${label}/
    mkdir -p rtt_plots_${label} && mv rtt_plots/* rtt_plots_${label}/
    mv log.txt log.${label}.txt
    """
}

process keep_chromosome {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/keep_chromosome"

    input:
    path assembly_dir

    output:
    path "${assembly_dir}.fasta"

    script:
    """
    Rscript $projectDir/bin/keep_chromosome.R $assembly_dir
    """
}

process remove_duplicates {
    container "$hgttree_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/remove_duplicates"
    
    input:
    path chromosomes
    
    output:
    path "chromosomes.nodup.fasta", emit: chromosomes_nodup
    path "duplicates.txt", emit: duplicates
    
    script:
    """
    touch duplicates.txt
    seqkit rmdup -s $chromosomes -D duplicates.txt -o chromosomes.nodup.fasta -j ${task.cpus}
    """
}

process root_tree {
	container "$r_container"
	storeDir "$launchDir/$params.resdir/root_tree_${label}"

	input:
	tuple val(label), path(snps), path(tree)

	output:
	tuple val(label), path(snps), path("rooted_trees.rds"), emit: rooted_trees
	path "rtt_metrics.${label}.rds"
	path "rtt_plots.${label}.pdf"
	path "log.${label}.txt"

	script:
	"""
	Rscript $projectDir/bin/root_tree.R \
		--project_dir $projectDir \
		--tree $tree \
		--assemblies $params.assemblies \
		--root_method $params.root_method \
		--root_topn $params.root_topn \
		--threads ${task.cpus}
	
	mv rooted_trees.rds rooted_trees.${label}.rds
	mv rtt_metrics.rds rtt_metrics.${label}.rds
	mv rtt_plots.pdf rtt_plots.${label}.pdf
	mv log.txt log.${label}.txt
	"""
}

process shrink_tree {
    container "$hgttree_container"
    storeDir "$launchDir/$params.resdir/shrink_tree"

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
    tuple path(snps), path("treeshrink.tre"), emit: shrinked_tree
    path snps, emit: snps
    path "treeshrink.txt"
    path "treeshrink_summary.txt"    

    script:
    """
    run_treeshrink.py --tree $tree --outprefix treeshrink --force --outdir .
    """
}


process prune_tree {
    container "$treepruner_container"
    storeDir "$launchDir/$params.resdir/treepruner_tree"

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
    tuple path(snps), path("treepruner.tree"), emit: pruned_tree
    path "treepruner.tree"

    script:
    """
    treepruner.py -i $tree -o treepruner.tree
    """
}

process snippy {
    container "$snippy_container"
    containerOptions "--bind $workDir:/scratch/tmp"
    storeDir "$launchDir/$params.resdir/snippy"

    input:
    tuple val(assembly_id), val(R1), val(R2), val(contigs), val(mode), val(error)
    each reference_genome

    output:
    path assembly_id

    script:
    if (mode == "paired")
        """
        snippy \
        --outdir $assembly_id \
        --ref $reference_genome \
        --R1 $R1 \
        --R2 $R2 \
        --force
        """
    else if (mode == "single")
        """
        snippy \
        --outdir $assembly_id \
        --ref $reference_genome \
        --se $R1 \
        --force
        """
    else if (mode == "contigs")
        """
        snippy \
        --outdir $assembly_id \
        --ref $reference_genome \
        --ctgs $contigs \
        --force
        """
    else
        error "Invalid mode: ${mode}"
}

process tidy_bootstrap_tree {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/tidy_bootstrap_tree"

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
    storeDir "$launchDir/$params.resdir/validate_input"

    output:
    path "assemblies.tsv"

    script:
    """
    Rscript $projectDir/bin/validate_input.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies
    """
}

workflow {
    // Validate input and create a list of genomes to process
    validate_input() | prep_snippy_input
    
    // Choose a reference genome
    if (params.reference_genome != null) {
        reference_genome = "$launchDir/$params.reference_genome"
        println("Reference genome for snippy: $reference_genome")
    } else {
        reference_genome = choose_reference_genome(validate_input.out)
    }
    
    // Construct pseudo-whole genomes and keep only chromosomes
    snippy_channel = prep_snippy_input.out.snippy_input | splitCsv(header: true, sep: "\t")

    snippy(snippy_channel, reference_genome) | keep_chromosome

    chromosomes = keep_chromosome.out.collectFile(
        name: "chromosomes.fasta",
        storeDir: "$launchDir/$params.resdir/"
    )

    // Remove duplicates, mask recombination, build tree, shrink, bootstrap

    chromosomes | remove_duplicates
    build_tree_ch = remove_duplicates.out.chromosomes_nodup | build_tree
    
    // Split the channel so both processes get the same input
    
    // Pruning
    shrink_tree_out = build_tree_ch | shrink_tree
    prune_tree_out  = build_tree_ch  | prune_tree

    // Root shrinked tree using mad, midpoint, rtt
    
    root_input = shrink_tree_out.shrinked_tree
    .combine(prune_tree_out.pruned_tree) { a, b ->
        tuple(a[0], a[1], b[1])
    }
    
    // Root both trees separately
    shrink_root_out = shrink_tree_out.shrinked_tree | root_tree
    prune_root_out  = prune_tree_out.pruned_tree   | root_tree


	// Root both trees separately
	shrink_root_out = shrink_tree_out.shrinked_tree.map { snps, tree -> tuple("shrink", snps, tree) } | root_tree
	prune_root_out  = prune_tree_out.pruned_tree.map { snps, tree -> tuple("prune", snps, tree) } | root_tree

	// Date both rooted trees separately
	shrink_date_out = shrink_root_out.rooted_trees | date_tree
	prune_date_out  = prune_root_out.rooted_trees | date_tree
	
	merged_dated_trees = date_tree_shrink.out.dated_trees.combine(
		date_tree_prune.out.dated_trees
	) { shrink_path, prune_path ->
		tuple(shrink_path, prune_path)
	}

    choose_dated_tree(merged_dated_trees)

    // Add tips that were removed as duplicates to final dated tree
    add_duplicates(
        choose_dated_tree.out.dated_big_tree,
        remove_duplicates.out.duplicates
    )
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Containers
gubbins_container = "stitam/prophyl:0.11"
hgttree_container = "mesti90/hgttree:2.11"
iqtree_container = "staphb/iqtree"
fasttree_container = "staphb/fasttree:latest"
r_container = "stitam/prophyl:0.11"
snippy_container = "staphb/snippy"

// Input parameters

params.assemblies = "$launchDir/assemblies.tsv"
params.genome_dir = "$launchDir/genomes"
params.reference_genome = null
params.gubbins_iterations = 10
params.bootstrap_replicates = 1000

// Output parameters

params.resdir = "results"

process add_duplicates {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/add_duplicates"

    input:
    path tree
    path duplicates

    output:
    path "dated_tree.rds"

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
    path dated_trees

    output:
    path "final_dated_tree.rds", emit: dated_big_tree
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/choose_dated_tree.R --trees $dated_trees
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
    --assemblies $assemblies \
    --genome_dir $params.genome_dir
    """
}

process create_genome_list {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/create_genome_list"

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
    container "$r_container"
    storeDir "$launchDir/$params.resdir/date_tree"

    input:
    tuple path(snps), path(rooted_trees)

    output:
    path "dated_trees.rds", emit: dated_trees
    path "rtt_plots/*.pdf"
    path "dated_trees/*.tre"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/date_tree.R \
    --project_dir $projectDir \
    --trees $rooted_trees \
    --snps $snps \
    --assemblies $params.assemblies \
    --threads ${task.cpus} \
    --branch_dimension snp_per_genome \
    --reroot false
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
    storeDir "$launchDir/$params.resdir/root_tree"

    input:
    tuple path(snps), path(shrinked_tree)

    output:
    tuple path(snps), path("rooted_trees.rds"), emit: rooted_trees
    // path "rooted_trees/*.tre"
    path "rtt_metrics.rds"
    path "rtt_plots.pdf"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/root_tree.R \
    --project_dir $projectDir \
    --tree $shrinked_tree \
    --assemblies $params.assemblies \
    --threads ${task.cpus}
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

process snippy_contig {
    container "$snippy_container"
    storeDir "$launchDir/$params.resdir/snippy_contig"

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
    storeDir "$launchDir/$params.resdir/snippy_paired"

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
    storeDir "$launchDir/$params.resdir/snippy_single"

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
        storeDir: "$launchDir/$params.resdir/"
    )

    // Remove duplicates, mask recombination, build tree, shrink, bootstrap

    chromosomes | remove_duplicates
    remove_duplicates.out.chromosomes_nodup | build_tree | shrink_tree //| bootstrap_tree | tidy_bootstrap_tree

    // Root shrinked tree using mad, midpoint, rtt
    shrink_tree.out.shrinked_tree | root_tree

    // Date all rooted trees
    date_tree(root_tree.out.rooted_trees)

    date_tree.out.dated_trees | choose_dated_tree

    // Add tips that were removed as duplicates to final dated tree
    add_duplicates(
        choose_dated_tree.out.dated_big_tree,
        remove_duplicates.out.duplicates
    )
}
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

as_targets = Channel.of("country", "k_serotype")
pastml_ch = Channel.fromPath("$launchDir/output/pastml/combined_ancestral_states.tab", checkIfExists: true)
tree_ch = Channel.fromPath("$launchDir/treedater_tree_with_time.nwk", checkIfExists: true)
meta_rda_ch = Channel.fromPath("$launchDir/treemeta.rda", checkIfExists: true)
meta_tsv_ch = Channel.fromPath("$launchDir/treemeta.tsv", checkIfExists: true)

process create_genome_list {
    container "stitam/r-packages"
    
    input:
    file "assemblies.tsv"

    output:
    file "log.txt"
    file "paired_reads.csv"
    file "single_reads.csv"
    file "contigs.csv"

    script:
    """
    Rscript $projectDir/bin/prep_snippy_input.R $params.genomedir
    """
}

process snippy_paired {
    container "staphb/snippy"

    input:
    tuple val(assembly_id), val(R1), val(R2)

    output:
    
    script:
    """
    snippy \
    --outdir "./jobs/${jobname}/output/snippy/paired/$R1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --R1 $R1 \
    --R2 $R2 \
    --force
    """
}

process predict_ancestral_states {
    container "evolbioinfo/pastml"
    
    input:
    tuple path(dated_tree), path(treemeta), val(target)

    output:
    tuple val(target), path("*")

    script:
    """
    pastml \
    -t $dated_tree \
    -d $treemeta \
    -c $target \
    --threads ${task.cpus} \
    --offline  \
    """
}

process simplify_ancestral_states {
    container "stitam/r-packages"

    input:
    tuple val(target), path(combined)

    output:
    path "marginals.rds"

    script:
    """
    Rscript $projectDir/bin/simplify_marginals.R $combined $target
    """
}

process prep_tree_tbl {
    container "stitam/r-packages"

    input:
    path tree
    path treemeta
    path marginals

    output:
    path "tree_tbl.rds"
    
    script:
    """
    Rscript $projectDir/bin/prep_tree_tbl.R $tree $treemeta $marginals
    """
}

process plot_tree {
    container "stitam/rplots:1.0"

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

workflow {
    tree_ch.combine(meta_tsv_ch).combine(as_targets) | predict_ancestral_states | simplify_ancestral_states
    // prep_tree_tbl(tree_ch, meta_ch, simplify_marginals.out) | plot_tree
}

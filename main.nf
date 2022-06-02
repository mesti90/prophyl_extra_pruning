#!/usr/bin/env nextflow

nextflow.enable.dsl=2

pastml_ch = Channel.fromPath("$launchDir/output/pastml/combined_ancestral_states.tab", checkIfExists: true)
tree_ch = Channel.fromPath("$launchDir/treedater_tree_with_time.nwk", checkIfExists: true)
meta_ch = Channel.fromPath("$launchDir/treemeta.rda", checkIfExists: true)

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

process simplify_marginals {
    container "stitam/r-packages"
    storeDir "$launchDir"

    input:
    path combined

    output:
    path "simplified_marginals.rds"

    script:
    """
    Rscript $projectDir/bin/simplify_marginals.R $combined
    """
}

process plot_tree {
    container "stitam/r-packages"

    input:
    path tree
    path treemeta
    path marginals

    output:
    path "dated_tree.pdf"
    path "dated_tree.png"
    
    script:
    """
    Rscript $projectDir/bin/plot_tree.R $tree $treemeta $marginals
    """
}

workflow {
    pastml_ch | simplify_marginals
    plot_tree(tree_ch, meta_ch, simplify_marginals.out)
}

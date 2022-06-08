#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

process predict_country {
    container "evolbioinfo/pastml"
    
    input:
    path dated_tree
    path treemeta

    output:
    path "*"

    script:
    """
    pastml \
    -t $dated_tree \
    -d $treemeta \
    -c "country" \
    --threads ${task.cpus} \
    --offline  \
    """
}

process simplify_country {
    container "stitam/r-packages"

    input:
    path combined

    output:
    path "country_marginals.rds"

    script:
    """
    Rscript $projectDir/bin/simplify_marginals.R $combined "country"
    """
}

process predict_k_serotype {
    container "evolbioinfo/pastml"
    
    input:
    path dated_tree
    path treemeta

    output:
    path "*"

    script:
    """
    pastml \
    -t $dated_tree \
    -d $treemeta \
    -c "k_serotype" \
    --threads ${task.cpus} \
    --offline  \
    """
}

process simplify_k_serotype {
    container "stitam/r-packages"

    input:
    path combined

    output:
    path "k_serotype_marginals.rds"

    script:
    """
    Rscript $projectDir/bin/simplify_marginals.R $combined "k_serotype"
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
    predict_country(tree_ch, meta_tsv_ch)
    predict_country.out | simplify_country
    predict_k_serotype(tree_ch, meta_tsv_ch)
    predict_k_serotype.out | simplify_k_serotype
    // prep_tree_tbl(tree_ch, meta_ch, simplify_marginals.out) | plot_tree
}

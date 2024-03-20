#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// CONTAINERS
pastml_container = "evolbioinfo/pastml"
r_container = "stitam/prophyl:0.11"

// PARAMETERS

// Output parameters

params.assemblies = "$launchDir/assemblies.tsv"
params.tree = "$launchDir/dated_tree.nwk"
params.target = "city"
params.method = "MPPA"
params.model = "JC"
params.resdir = "results"

// Input files

// Processes 

process predict_ancestral_states {
    container "$pastml_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/predict_ancestral_states"

    input:
    val(target)

    output:
    tuple val(target), path(target)

    script:
    """
    pastml \
    -t $params.tree \
    -d $params.assemblies \
    -c $target \
    --prediction_method $params.method \
    -m $params.model \
    --threads ${task.cpus} \
    --work_dir $params.target \
    --html_compressed tree_with_ancestral_states.html \
    --offline
    """
}

process prep_tree_tbl {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/prep_tree_tbl"

    input:
    path ancestral_states

    output:
    path "tree_tbl.rds"
    path "tree_tbl.tsv"
    
    script:
    """
    Rscript $projectDir/bin/prep_tree_tbl.R \
    --project_dir $projectDir \
    --assemblies $params.assemblies \
    --tree $params.tree \
    --ancestral_states $ancestral_states
    """
}

process tidy_ancestral_states {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/tidy_ancestral_states"

    input:
    tuple val(target), path(combined)

    output:
    path "${target}.tsv"

    script:
    """
    Rscript $projectDir/bin/tidy_ancestral_states.R $combined $target
    """
}

// Workflow

workflow {
    // predict ancestral states for the target variable
    predict_ancestral_states(params.target)
    // convert ancestral states to a table
    predict_ancestral_states.out | tidy_ancestral_states
    // combine with other metadata count state changes etc.
    tidy_ancestral_states.out | prep_tree_tbl
}
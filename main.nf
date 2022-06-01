#!/usr/bin/env nextflow

params.genomedir="/home/tamas/Server/hpc/node8_R10/kintses_lab/aci"

nextflow.enable.dsl=2

assembly_ch = Channel.fromPath("$launchDir/assemblies.tsv", checkIfExists: true)

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

workflow {
    create_genome_list(assembly_ch)
    create_genome_list.out[1].splitCsv().view()
//    format_mlst(predict_mlst(assembly_tuple_ch))
//    predict_k_serotype(assembly_tuple_ch)
//    predict_o_serotype(assembly_tuple_ch) 
//    format_kaptive(predict_k_serotype.out.join(predict_o_serotype.out))
//    format_resgenes(predict_resgenes(predict_orfs(assembly_tuple_ch)))  
//    bind_predictions(format_mlst.out.join(format_kaptive.out.join(format_resgenes.out)))

//    bind_predictions.out.collectFile(name: "sample.tsv", newLine: false, keepHeader: true, storeDir: "$launchDir")
   }

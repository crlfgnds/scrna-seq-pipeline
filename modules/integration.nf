#!/usr/bin/env nextflow

//copy inputs, outputs and script from 05_integration.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process integration {

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path rds_file
    val integration_features
    val method

    output:
    path 'seurat_integrated.rds'

    script:
    """
    Rscript ${projectDir}/scripts/05_integration.R ${rds_file} ${integration_features} ${method}

    """
}

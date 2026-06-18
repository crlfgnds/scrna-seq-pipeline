#!/usr/bin/env nextflow

//copy inputs, outputs and script from 04_normalize_hvg.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process normalize_hvg {

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path rds_file
    val n_features
    val split_by

    output:
    path 'seurat_list_normalized.rds'
    path 'hvg_plot.pdf'

    script:
    """
    Rscript ${projectDir}/scripts/04_normalize_hvg.R ${rds_file} ${n_features} ${split_by}

    """
}

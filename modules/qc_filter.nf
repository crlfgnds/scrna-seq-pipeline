#!/usr/bin/env nextflow

//copy inputs, outputs and script from 02_qc_filter.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process qc_filter {

    publishDir "${params.outdir}/preprocessing/qc", mode: 'copy'

    input:
    path rds_file
    val mito_thresh
    val min_features
    val max_features
    val min_counts
    val max_counts

    output:
    path 'seurat_filtered.rds'
    path 'qc_before_filtering.pdf'
    path 'qc_after_filtering.pdf'

    script:
    """
    Rscript ${projectDir}/scripts/02_qc_filter.R ${rds_file} ${mito_thresh} ${min_features} ${max_features} ${min_counts} ${max_counts}

    """
}

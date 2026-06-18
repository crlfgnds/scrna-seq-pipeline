#!/usr/bin/env nextflow

process hto_demux {

    publishDir "${params.outdir}/preprocessing/hto", mode: 'copy'

    input:
    path rds_file
    path hto_h5
    val hashtags
    val hashtag_labels
    val group_pattern
    val ident1
    val ident2
    val experiment_col

    output:
    path 'seurat_demuxed.rds'
    path 'hto_ridge.pdf'

    script:
    """
    Rscript ${projectDir}/scripts/03_hto_demux.R ${rds_file} ${hto_h5} ${hashtags} ${hashtag_labels} ${group_pattern} ${ident1} ${ident2} ${experiment_col}
    """
}

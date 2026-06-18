#!/usr/bin/env nextflow

//copy inputs, outputs and script from 01_create_seurat.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process create_seurat {

    publishDir "${params.outdir}/preprocessing", mode: 'copy'

    input:
    path h5_file
    val sample_name
  
    output:
    path 'seurat_raw.rds'

    script:
    """
    Rscript ${projectDir}/scripts/01_create_seurat.R ${h5_file} ${sample_name}

    """
}

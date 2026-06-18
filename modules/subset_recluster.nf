#!/usr/bin/env nextflow

//copy inputs, outputs and script from 08_subset_recluster.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process subset_recluster {

    publishDir "${params.outdir}/subset", mode: 'copy'

    input:
    path rds_file
    val clusters
    val subset_name
    val experiment_col
    val experiments
    val subset_n_dims
    val subset_resolution
    val subset_palette_arg

    output:
    path "${subset_name}_subsetted.rds"
    path "${subset_name}_elbow.pdf"
    path "${subset_name}_umap_resolutions.pdf"
    path "${subset_name}_umap_final.pdf"

    script:
    """
  Rscript ${projectDir}/scripts/08_subset_recluster.R ${rds_file} ${clusters} ${subset_name} ${experiment_col} ${experiments} ${subset_n_dims} ${subset_resolution} ${subset_palette_arg}
    """
}

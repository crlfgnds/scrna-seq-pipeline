#!/usr/bin/env nextflow

//copy inputs, outputs and script from 06_cluster.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process cluster {

    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path rds_file
    val n_dims
    val resolution
    val n_neighbors
    val min_dist
    val group_by_var
    val palette_arg

    output:
    path 'seurat_clustered.rds'
    path 'elbow_plot.pdf'
    path 'umap_resolutions.pdf'
    path 'umap_clusters.pdf'
    path 'umap_by_group.pdf'

    script:
    """
    Rscript ${projectDir}/scripts/06_cluster.R ${rds_file} ${n_dims} ${resolution} ${n_neighbors} ${min_dist} ${group_by_var} ${palette_arg}

    """
}

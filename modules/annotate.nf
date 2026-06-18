#!/usr/bin/env nextflow

//copy inputs, outputs and script from 07_annotate.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process annotate {

    publishDir "${params.outdir}/annotation", mode: 'copy'

    input:
    path rds_file
    path marker_file
    path annotation_file
    val min_pct
    val logfc_thresh
    val top_n_genes
    val palette_arg

    output:
    path 'seurat_annotated.rds'
    path 'all_markers.xlsx'
    path 'top_genes_per_cluster.xlsx'
    path 'heatmap_top_markers.pdf'
    path 'umap_annotated.pdf'
    path 'dotplot_markers.pdf',     optional: true
    path 'featureplot_markers.pdf', optional: true
    path 'violin_markers.pdf',      optional: true

    script:
    """
    Rscript ${projectDir}/scripts/07_annotate.R ${rds_file} ${marker_file} ${annotation_file} ${min_pct} ${logfc_thresh} ${top_n_genes} ${palette_arg} 

    """
}

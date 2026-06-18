#!/usr/bin/env nextflow

process cellchat {

    publishDir "${params.outdir}/downstream/cellchat", mode: 'copy'

    input:
    path rds_file
    val group_col
    val group_levels
    val celltype_col
    val organism
    val source_cells

    output:
    path "cellchat_compare_interactions.pdf"
    path "cellchat_diffInteraction_circles.pdf"
    path "cellchat_heatmap.pdf"
    path "cellchat_circle_split.pdf"
    path "cellchat_circle_per_celltype.pdf"
    path "cellchat_signalingRole_scatter.pdf"
    path "cellchat_signalingRole_heatmap_incoming.pdf"
    path "cellchat_signalingRole_heatmap_outgoing.pdf"
    path "cellchat_rankNet_*.pdf"
    path "cellchat_bubble_*.pdf"
    path "cellchat_communication_*.xlsx"
    path "cellchat_*.rds"

    script:
    """
    Rscript ${projectDir}/scripts/13_cellchat.R ${rds_file} ${group_col} ${group_levels} ${celltype_col} ${organism} ${source_cells}
    """
}

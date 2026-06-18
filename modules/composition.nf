#!/usr/bin/env nextflow

process composition {

    publishDir "${params.outdir}/downstream/composition", mode: 'copy'

    input:
    path rds_file
    val celltype_col
    val group_col
    val sample_col
    val ident1
    val ident2
    val palette_arg
    val time_col
    val time_levels

    output:
    path "composition_boxplot.pdf"
    path "composition_barplot.pdf"
    path "composition_per_sample.xlsx"
    path "composition_alluvial.pdf", optional: true

    script:
    """
    Rscript ${projectDir}/scripts/11_composition.R ${rds_file} ${celltype_col} ${group_col} ${sample_col} ${ident1} ${ident2} ${palette_arg} ${time_col} ${time_levels}
    """
}

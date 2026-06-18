#!/usr/bin/env nextflow

//copy inputs, outputs and script from 09_de.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process de {

    publishDir "${params.outdir}/downstream/de", mode: 'copy'

    input:
    path rds_file
    val group_col
    val ident1
    val ident2
    val min_pct
    val logfc_thresh
    val padj_cutoff
    val volcano_xlim
    val volcano_fc
    val volcano_p

    output:
    path 'DE_whole_object.xlsx'
    path 'DE_all_results.xlsx'
    path 'DE_significant_results.xlsx'
    path "volcano_*.pdf"

    script:
    """
  Rscript ${projectDir}/scripts/09_de.R ${rds_file} ${group_col} ${ident1} ${ident2} ${min_pct} ${logfc_thresh} ${padj_cutoff} ${volcano_xlim} ${volcano_fc} ${volcano_p}
    """
}

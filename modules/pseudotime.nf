#!/usr/bin/env nextflow

//copy inputs, outputs and script from 12_pseudotime.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process pseudotime {

    publishDir "${params.outdir}/downstream/pseudotime", mode: 'copy'

    input:
    path rds_file
    val root_cluster
    val group_col
    val group_levels
    val palette_arg

    output:
    path "pseudotime_all.pdf"
    path "pseudotime_boxplot.pdf"
    path "pseudotime_split_groups.pdf"
    path "cds_all.rds"
    path "cds_list_groups.rds"

    script:
    """
  Rscript ${projectDir}/scripts/12_pseudotime.R ${rds_file} ${root_cluster} ${group_col} ${group_levels} ${palette_arg} 
    """
}

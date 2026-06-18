#!/usr/bin/env nextflow

//copy inputs, outputs and script from 10_go.R
//inputs are the args
//output is the rds file name
//script is Rscript + name of script + args

process go {

    publishDir "${params.outdir}/downstream/go", mode: 'copy'

    input:
    path de_file
    path whole_de_file
    val organism
    val ont
    val n_top_terms
    val min_gs_size
    val max_gs_size
    val score_type
    val pvalue_cutoff
    val direction 

    output:
    path "GO_*_barplot.pdf"
    path "GO_*.xlsx"

    script:
    """
  Rscript ${projectDir}/scripts/10_go.R ${de_file} ${whole_de_file} ${organism} ${ont} ${n_top_terms} ${min_gs_size} ${max_gs_size} ${score_type} ${pvalue_cutoff} ${direction} 
    """
}

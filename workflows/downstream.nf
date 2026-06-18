#!/usr/bin/env nextflow

include { de }          from '../modules/de.nf'
include { go }          from '../modules/go.nf'
include { composition } from '../modules/composition.nf'
include { pseudotime }  from '../modules/pseudotime.nf'
include { cellchat }    from '../modules/cellchat.nf'

workflow DOWNSTREAM {
    take:
    rds_file
    group_col
    ident1
    ident2
    min_pct
    logfc_thresh
    padj_cutoff
    volcano_xlim
    volcano_fc
    volcano_p
    organism
    ont
    n_top_terms
    min_gs_size
    max_gs_size
    score_type
    pvalue_cutoff
    direction
    celltype_col
    sample_col
    palette_arg
    time_col
    time_levels
    root_cluster
    group_levels
    source_cells
    run_de
    run_composition
    run_pseudotime
    run_cellchat

    main:
    // de → go are a unit: go consumes de's output table
    if (run_de) {
        de(rds_file, group_col, ident1, ident2, min_pct, logfc_thresh, padj_cutoff, volcano_xlim, volcano_fc, volcano_p)
        go(de.out[1], de.out[0], organism, ont, n_top_terms, min_gs_size, max_gs_size, score_type, pvalue_cutoff, direction)
    }
    if (run_composition) composition(rds_file, celltype_col, group_col, sample_col, ident1, ident2, palette_arg, time_col, time_levels)
    if (run_pseudotime)  pseudotime(rds_file, root_cluster, group_col, group_levels, palette_arg)
    if (run_cellchat)    cellchat(rds_file, group_col, group_levels, celltype_col, organism, source_cells)
}

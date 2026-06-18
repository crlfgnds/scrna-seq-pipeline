#!/usr/bin/env nextflow

include { PREPROCESSING } from './workflows/preprocessing.nf'
include { cluster }          from './modules/cluster.nf'
include { annotate }         from './modules/annotate.nf'
include { subset_recluster } from './modules/subset_recluster.nf'
include { DOWNSTREAM }    from './workflows/downstream.nf'

//subworkflow (uppercase). Finishes and saves to disk. 
workflow PHASE1 {
    PREPROCESSING(
        params.h5_file,
        params.sample_name,
        params.mito_thresh,
        params.min_features,
        params.max_features,
        params.min_counts,
        params.max_counts,
        params.hto_h5,
        params.hashtags,
        params.hashtag_labels,
        params.group_pattern,
        params.ident1,
        params.ident2,
        params.experiment_col,
        params.n_features,
        params.split_by,
        params.integration_features,
        params.integration_method,
        params.hto
    )
    cluster(PREPROCESSING.out.rds,
            params.n_dims,
            params.resolution,
            params.n_neighbors,
            params.min_dist,
            params.group_by_var,
            params.palette_arg)
}
// CLUSTER_ONLY: enter here with an already-integrated object (skip preprocessing).
// Reads rds from disk — point --rds_file at your integrated Seurat object.
workflow CLUSTER_ONLY {
    rds_ch = channel.value(file(params.rds_file))
    cluster(rds_ch,
            params.n_dims,
            params.resolution,
            params.n_neighbors,
            params.min_dist,
            params.group_by_var,
            params.palette_arg)
}
// PHASE2 calls the annotate module (lowercase). Separate run → reads rds from disk.
workflow PHASE2{
    rds_ch = channel.value(file(params.rds_file))
    annotate(
        rds_ch,
        params.marker_file,
        params.annotation_file,
        params.min_pct,
        params.logfc_thresh,
        params.top_n_genes,
        params.palette_arg
    )
}
// PHASE3 calls the subset_recluster module (lowercase). Optional. Reads rds from disk.
workflow PHASE3{
    rds_ch = channel.value(file(params.rds_file))
    subset_recluster(
        rds_ch,
        params.clusters,
        params.subset_name,
        params.experiment_col,
        params.experiments,
        params.subset_n_dims,
        params.subset_resolution,
        params.subset_palette_arg
    )
}
// PHASE4 calls the DOWNSTREAM sub-workflow (uppercase). Reads rds from disk —
// point --rds_file at the annotated OR the subsetted object.
workflow PHASE4{
    rds_ch = channel.value(file(params.rds_file))
    DOWNSTREAM(
        rds_ch,
        params.group_col,
        params.ident1,
        params.ident2,
        params.min_pct,
        params.logfc_thresh,
        params.padj_cutoff,
        params.volcano_xlim,
        params.volcano_fc,
        params.volcano_p,
        params.organism,
        params.ont,
        params.n_top_terms,
        params.min_gs_size,
        params.max_gs_size,
        params.score_type,
        params.pvalue_cutoff,
        params.direction,
        params.celltype_col,
        params.sample_col,
        params.palette_arg,
        params.time_col,
        params.time_levels,
        params.root_cluster,
        params.group_levels,
        params.source_cells,
        params.run_de,
        params.run_composition,
        params.run_pseudotime,
        params.run_cellchat
    )
}

// ENTRY POINT — one workflow, selects which phase to run via --step.
// Nextflow 26's strict parser does NOT support -entry, so phases are dispatched by parameter:
//   nextflow run main.nf --step phase1 -profile local -params-file params/x.yaml
workflow {
    if      ( params.step == 'phase1'       ) PHASE1()
    else if ( params.step == 'cluster_only' ) CLUSTER_ONLY()
    else if ( params.step == 'phase2'       ) PHASE2()
    else if ( params.step == 'phase3'       ) PHASE3()
    else if ( params.step == 'phase4'       ) PHASE4()
    else error "Unknown --step '${params.step}'. Use one of: phase1 | cluster_only | phase2 | phase3 | phase4"
}

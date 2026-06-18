#!/usr/bin/env nextflow

include { create_seurat }  from '../modules/create_seurat.nf'
include { qc_filter }      from '../modules/qc_filter.nf'
include { hto_demux }      from '../modules/hto_demux.nf'
include { normalize_hvg }  from '../modules/normalize_hvg.nf'
include { integration }    from '../modules/integration.nf'

workflow PREPROCESSING {
    take:
    h5_file
    sample_name
    mito_thresh
    min_features
    max_features
    min_counts
    max_counts
    hto_h5
    hashtags
    hashtag_labels
    group_pattern
    ident1
    ident2
    experiment_col
    n_features
    split_by
    integration_features
    integration_method
    hto

    main:
    // .out[0] always = the .rds; later indices are plots. Index explicitly so adding
    // plot outputs never silently re-wires the chain.
    create_seurat(h5_file, sample_name)
    qc_filter(create_seurat.out[0], mito_thresh, min_features, max_features, min_counts, max_counts)

    if (hto) {
        hto_demux(qc_filter.out[0], hto_h5, hashtags, hashtag_labels, group_pattern, ident1, ident2, experiment_col)
        normalize_hvg(hto_demux.out[0], n_features, split_by)
    } else {
        normalize_hvg(qc_filter.out[0], n_features, split_by)
    }

    integration(normalize_hvg.out[0], integration_features, integration_method)

    emit:
    rds = integration.out[0]
}

# scRNA-seq Analysis Pipeline

A modular, reproducible single-cell RNA-seq analysis pipeline built with **Nextflow (DSL2)** and **R/Seurat**, fully containerised with Docker/Singularity.

It takes raw 10x Genomics output through QC, hashtag demultiplexing, normalisation, integration, clustering, annotation, differential expression, GO enrichment, composition analysis, pseudotime, and cell–cell communication (CellChat) — with explicit, documented decision points where a human must inspect results before continuing.

---

## Design philosophy

scRNA-seq analysis is **not** a press-play-and-walk-away process. Several steps require you to look at a plot and make a judgement call:

- How many principal components to use for clustering (elbow plot)
- What clustering resolution to choose
- Which cell-type names to assign to each cluster
- Whether to subset and re-cluster a population
- Where to set volcano-plot thresholds

This pipeline makes those decision points **first-class**. Instead of hiding the judgement calls, it splits the workflow into four **phases**, each ending exactly where you need to stop, inspect output, and feed a decision back in. Every decision becomes a parameter (a YAML value, a CSV file, or an `.rds` path) — so the analysis stays reproducible and re-runnable while you remain in control of the science.

```
PHASE1  preprocess + cluster   →  stop: choose n_dims, resolution
PHASE2  annotate               →  stop: assign cell-type names
PHASE3  subset_recluster (opt) →  stop: decide whether to subset
PHASE4  downstream analysis    →  DE, GO, composition, pseudotime, CellChat
```

---

## Architecture

Three layers, cleanly separated:

```
main.nf                    entry point — phases selected via --step
  │
  ├── workflows/           sub-workflows = chains of modules (UPPERCASE)
  │     PREPROCESSING, DOWNSTREAM
  │
  └── modules/             one process each = one R script (lowercase)
          │
          └── scripts/     the actual R analysis code (01..13)
```

- **`modules/*.nf`** — each wraps a single R script, declaring its inputs and output files.
- **`workflows/*.nf`** — chain modules together, passing outputs to inputs.
- **`main.nf`** — defines the four phase entry points and feeds them parameters.
- **`nextflow.config`** — all parameters with defaults, profiles, container.
- **`conf/`** — profile-specific settings (`local.config` = Docker, `hpc.config` = Singularity + SLURM).
- **`params/`** — dataset-specific parameter files (override the defaults).

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.10
- **Docker** (local) or **Singularity/Apptainer** (HPC)

All R dependencies (Seurat, CellChat, monocle3, clusterProfiler, EnhancedVolcano, etc.) are baked into the container image — nothing to install locally.

**Container:** [`crlfgnds/scrna-seq-pipeline:latest`](https://hub.docker.com/r/crlfgnds/scrna-seq-pipeline) (built on `bioconductor/bioconductor_docker:RELEASE_3_19`).

---

## Usage

Each phase is its own command. You select a phase with `--step`, a compute profile with `-profile`, and a dataset with `-params-file`.

> Nextflow 26's strict parser does not support `-entry`; phases are selected with the `--step` parameter instead.

### Phase 1 — preprocess + cluster
```bash
nextflow run main.nf --step phase1 \
  -profile local \
  -params-file params/neutrophil_mouse.yaml
```
Then inspect `results/clustering/elbow_plot.pdf` and `umap_resolutions.pdf`. Set `n_dims` and `resolution` in the YAML.
Re-run with `-resume` — only the `cluster` step recomputes, everything upstream is cached.

### Phase 2 — annotate
```bash
nextflow run main.nf --step phase2 \
  -profile local \
  -params-file params/neutrophil_mouse.yaml \
  --rds_file results/clustering/seurat_clustered.rds
```
Inspect the marker tables/plots, write a `cluster,cell_type` CSV, set it as `annotation_file` in the YAML, re-run.

### Phase 3 — subset & re-cluster (optional)
```bash
nextflow run main.nf --step phase3 \
  -profile local \
  -params-file params/neutrophil_mouse.yaml \
  --rds_file results/annotation/seurat_annotated.rds
```

### Phase 4 — downstream analysis
```bash
nextflow run main.nf --step phase4 \
  -profile local \
  -params-file params/neutrophil_mouse.yaml \
  --rds_file results/annotation/seurat_annotated.rds   # or a results/subset/*_subsetted.rds
```
Point `--rds_file` at the **whole annotated object** or a **subsetted population** — whichever you want to analyse.

> **HPC:** swap `-profile local` for `-profile hpc` to run with Singularity on a SLURM cluster (adjust queue/resources in `conf/hpc.config`).

### Starting from an existing Seurat object

Because every phase after PHASE1 reads its input from `--rds_file`, you can enter the pipeline wherever your object already is — no need to start from raw `.h5`:

| Object you already have | Enter at |
|---|---|
| Integrated, not yet clustered | `--step cluster_only` |
| Already clustered | `--step phase2` (annotate) |
| Already annotated | `--step phase3` or `--step phase4` |

```bash
# e.g. cluster an already-integrated object, skipping all preprocessing
nextflow run main.nf --step cluster_only -profile local \
  -params-file params/neutrophil_mouse.yaml \
  --rds_file my_integrated_object.rds
```

> **What the object must already contain depends on where you enter** — because jumping in late means you skipped the step that would have created it:
> - `--step cluster_only` → needs an `integrated` assay (clustering runs on it)
> - `--step phase2` → needs cluster identities (i.e. it's already been clustered)
> - `--step phase3` / `--step phase4` → needs a `cell_type` column and an RNA assay (normally created by `annotate` in phase2)
>
> In the full flow from `.h5`, you never supply these — the pipeline creates `cell_type` during annotation (phase2) and the assays during preprocessing.

### Skipping downstream analyses

PHASE4 runs DE+GO, composition, pseudotime and CellChat by default. Toggle any off in the YAML (or on the CLI):

```yaml
run_pseudotime: false
run_cellchat:   false
```

(`run_de` controls DE **and** GO, since GO consumes the DE table.)

---

## The inspect → tune → resume loop

For parameter decisions (PCs, resolution, volcano limits), you don't restart from scratch. Nextflow caches every step by a hash of its inputs:

1. Run the phase with default/guessed parameters.
2. Inspect the output plots.
3. Change the relevant value in the YAML.
4. Re-run with `-resume` — only the affected step (and anything downstream of it) recomputes.

The expensive steps (integration especially) run **once** and are reused.

---

## Pipeline steps & key outputs

| # | Module | Purpose | Key outputs |
|---|--------|---------|-------------|
| 01 | create_seurat | Build Seurat object from 10x `.h5` | `seurat_raw.rds` |
| 02 | qc_filter | Filter cells on mito/counts/features | `seurat_filtered.rds`, QC plots |
| 03 | hto_demux | Hashtag demultiplexing (optional) | `seurat_demuxed.rds`, ridge plot |
| 04 | normalize_hvg | Normalise + highly variable genes | `seurat_list_normalized.rds`, HVG plot |
| 05 | integration | Integrate samples (CCA/RPCA/Harmony) | `seurat_integrated.rds` |
| 06 | cluster | PCA, UMAP, clustering | `seurat_clustered.rds`, elbow + UMAP plots |
| 07 | annotate | Cluster markers + cell-type labels | `seurat_annotated.rds`, marker tables + plots |
| 08 | subset_recluster | Subset a population & re-cluster | `*_subsetted.rds`, UMAPs |
| 09 | de | Differential expression (per cluster + whole) | DE tables, volcano plots |
| 10 | go | GO enrichment (up & down) | GO bar plots + tables |
| 11 | composition | Cell-type frequencies between conditions | boxplot, barplot, alluvial |
| 12 | pseudotime | Trajectory inference (monocle3) | trajectory plots, CDS objects |
| 13 | cellchat | Cell–cell communication | CellChat figures + tables |

---

## Parameters

All parameters and their defaults live in `nextflow.config`. Dataset-specific values are set in a YAML under `params/` (see `params/neutrophil_mouse.yaml` and `params/human_skin.yaml` for worked examples). Any parameter can also be overridden on the command line, e.g. `--n_dims 12`.

---

## Author

**Anna Carolina Fagundes** — computational immunologist.
Developed from PhD single-cell analyses of intestinal immune populations.

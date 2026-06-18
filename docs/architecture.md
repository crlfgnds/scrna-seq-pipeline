# Architecture — how the pieces fit (and why it looks repetitive)

A personal reference for understanding *why* the pipeline is structured the way it is.
Nextflow is verbose, but almost nothing is truly redundant — each layer is a **border
between two different runtimes**, and a value has to be "stamped through customs" at each one.

---

## 1. The three runtimes

The pipeline spans three separate worlds that don't share variables:

```
   NEXTFLOW  ──►  SHELL (bash)  ──►  R
   (orchestration)  (the command)    (the analysis)
```

A parameter like `n_dims` starts life in Nextflow, gets injected into a shell command,
and is finally read by an R script. Each world names it in its own language — that's why
the same name appears several times. It's **translation**, not repetition.

---

## 2. One parameter, traced end to end

Follow `n_dims` (number of PCs for clustering) through every layer:

```
LAYER                              CODE                              ITS ONE JOB
────────────────────────────────  ────────────────────────────────  ──────────────────────────
nextflow.config                    params { n_dims = 17 }            define it + default value
                                                                     (single source of truth)

main.nf  (workflow PHASE1)         cluster(.., params.n_dims, ..)    hand the VALUE to the call

modules/cluster.nf  input:         val n_dims                        process declares it accepts
                                                                     a value here  (Nextflow)

modules/cluster.nf  script:        Rscript .. ${n_dims}              inject it into the shell
                                                                     command  (bash)

scripts/06_cluster.R               n_dims <- as.integer(args[2])     R reads it from the CLI  (R)
```

The bottom three *look* redundant but are three different languages handing the value across
two borders: **Nextflow → bash → R**. Remove any one and the chain breaks.

---

## 3. What's load-bearing vs what's a design choice

```
LOAD-BEARING (cannot be removed):
  • config definition            no default, no single source without it
  • module `input:` declaration  Nextflow must know the process signature
  • module `script: ${}`         the only way into the shell command
  • R `args[]`                   R's only way to receive CLI input

DELIBERATE DESIGN (adds layers, buys something):
  • sub-workflow `take:` blocks   makes PREPROCESSING / DOWNSTREAM reusable, self-contained
  • the phase split (--step)      enables the stop → inspect → resume workflow
```

If it ever feels like "too many steps," it's almost always one of the *design* layers —
and those exist to buy modularity and the human-in-the-loop decision points, not for ceremony.

---

## 4. The layer cake (who calls whom)

```
nextflow.config        all params + defaults + profiles + container
        │
main.nf                ENTRY POINT — one workflow{}, dispatches by --step
        │                 selects PHASE1 | CLUSTER_ONLY | PHASE2 | PHASE3 | PHASE4
        │
        ├── workflows/   SUB-WORKFLOWS (UPPERCASE) — chains of modules
        │      PREPROCESSING, DOWNSTREAM
        │              │
        └── modules/   one process each (lowercase) = wraps one R script
                   │      declares input:, output:, publishDir, script:
                   │
                   └── scripts/   the actual R analysis (01..13)
```

Naming rule: **UPPERCASE = a group of steps (sub-workflow). lowercase = one step (module/process).**

---

## 5. Why phases read from disk (`channel.value(file(params.rds_file))`)

Each phase is a **separate command run at a separate time** (so you can stop and make a
decision in between). When PHASE2 runs, PHASE1 is long over and its in-memory channels are
gone — the only thing connecting them is the `.rds` file saved on disk. So PHASE2 must pick
that file back up:

```
SAME run  (inside PHASE1):  cluster(PREPROCESSING.out.rds)   ← live channel, .out
NEW run   (PHASE2 later):   channel.value(file(params.rds_file))  ← reload from disk
```

`channel.value(...)` (not `channel.fromPath`) so the file can be read by **multiple**
processes in that phase without being consumed once and emptied.

---

## 6. Why `.out[0]` and not `.out`

When a process has **one** output, `proc.out` *is* that file. When it has **several**
(e.g. an `.rds` plus plots), `proc.out` becomes an indexed list — `.out[0]` is the rds,
`.out[1]` a plot, etc. Always index explicitly when chaining, so adding a plot output later
never silently re-wires the chain to the wrong file.

---

## 7. Where outputs land

Each module has a `publishDir`, so its outputs are copied into `results/`:

```
results/
├── preprocessing/        seurat objects, QC + HVG plots
│   ├── qc/               qc_before/after_filtering.pdf
│   └── hto/              demux ridge plot (if --hto)
├── clustering/           elbow, resolution + cluster UMAPs   ← inspect to pick n_dims
├── annotation/           marker tables, heatmap, feature/violin/dot plots
├── subset/               subsetted object + its UMAPs (optional)
└── downstream/
    ├── de/               DE tables + volcano plots
    ├── go/               GO bar plots + tables
    ├── composition/      boxplot, barplot, alluvial
    ├── pseudotime/       trajectory plots + CDS objects
    └── cellchat/         CellChat figures + tables
```

# ContextNet Tutorials

ContextNet extends NicheNet by replacing the generic gene regulatory network with a context-specific gene regulatory network (GRN). Use ContextNet when you have a regulatory network inferred from your own single-cell or multiomic system and want ligand-target predictions to reflect that context rather than a generic prior. The tutorials are organized around reader questions: what input is needed, what biological question is being asked, what object is produced, and how to interpret the output.

## Tutorial Map

| Tutorial | Biological or training question | Data mode | Main output | Use when |
|----------|---------------------------------|-----------|-------------|----------|
| [01 Getting Started](01-getting-started.ipynb) / [Rmd](01-getting-started.Rmd) | Does my installation work, and what do ContextNet inputs and outputs look like? | Packaged synthetic smoke-test fixture | `ContextGRN`, `ContextNetwork`, `ContextResult`, expected plots | You are learning the object model or checking an install |
| [02 Standard Workflow](02-standard-workflow.ipynb) / [Rmd](02-standard-workflow.Rmd) | Which sender ligands may explain basal, parabasal, suprabasal, and dividing receiver states within contralateral tissue? | Full local Seurat + context-specific GRN | receiver-specific ligand activity, target heatmaps, ccflowR Sankey river paths | You want a single-tissue ContextNet analysis across epithelial receiver states |
| [03 Multi-Condition Workflow](03-multicondition-workflow.ipynb) / [Rmd](03-multicondition-workflow.Rmd) | Which communication programs change reproducibly across samples? | Full local replicate-aware MultiNicheNet analysis | prioritization table, abundance information, ccflowR Sankey paths | You have biological replicates and need sample-aware inference |

RStudio users can open the companion `.Rmd` files in this directory. They contain the same source cells as the notebooks and are meant for interactive execution, not package build-time vignettes.

## Running the Tutorials

- **Installation.** Install ContextNet first using the repository `README.md` or the micromamba guide. The tutorials do not install packages; they start from biological or training questions and use setup cells only to load libraries and local data.

- **Kernel.** The notebooks use an R kernel named `contextnet_r`, displayed in Jupyter as **ContextNet R**. Register it once with `IRkernel::installspec(name = "contextnet_r", displayname = "ContextNet R")`. If Jupyter also shows an old raw `ir` kernel for the same environment, remove the confusing entry with `jupyter kernelspec remove ir -f` after confirming that `contextnet_r` works.

- **Working directory.** Run a notebook from the repository root or from the `notebooks/` directory. The setup cell normalizes the working directory to the repository root.

- **ccflowR river plots.** ccflowR is optional visualization support for Sankey/river plots. It helps visualize high-priority sender -> ligand -> receptor -> receiver paths after ContextNet or MultiNicheNet has ranked communication hypotheses. Install it into the same **ContextNet R** environment once the public `settylab/ccflowR` repository is available. Without ccflowR, core ContextNet analyses still run, but Sankey/river plot cells should be skipped.

- **Outputs and figures.** Tutorial 01 is committed with outputs from the packaged smoke-test fixture. Tutorials 02 and 03 are committed with current executed outputs from the local HNSCC reference analyses so readers can see the expected tables and figures. Re-run notebooks whenever inputs, package versions, or plotting code change.

## Data Requirements

Tutorial 01 uses only the packaged public fixture under `inst/extdata/contextnet_smoke/`. It is synthetic and exists only to test installation and object interpretation.

Tutorials 02 and 03 are full local-data tutorials. To run them as written you need, under `data/` or through the documented environment variables:

| Input | File used by the tutorials | Purpose |
|-------|----------------------------|---------|
| Context-specific GRN | `gene_x_tf.csv` | target-gene x transcription-factor regulatory weights from a GRN inference workflow |
| Seurat object | `seurat_obj.rds` | expression matrix plus cell-type, condition, and sample metadata |
| NicheNet L-R prior | `lr_network_human_21122021.rds` | ligand-receptor interactions |
| NicheNet weighted prior | `weighted_networks_nsga2r_final.rds` | signaling and regulatory priors for propagation |
| NicheNet L-T matrix (optional) | `ligand_target_matrix_nsga2r_final.rds` | reference only; not used by ContextNet integration |

The NicheNet priors are distributed with the official NicheNet v2 tutorials and can be fetched once with `load_nichenet_data(cache_dir = "data", version = "v2")`. The full HNSCC Seurat object is not bundled because it is large and study-specific.

## How to Interpret the Object Model

ContextNet separates the analysis into three object classes:

- `ContextGRN`: regulatory evidence from the profiled biological context. It stores transcription-factor to target-gene edges imported from a GRN inference workflow.
- `ContextNetwork`: reusable propagated signaling model. It stores the context-specific ligand-target matrix created by integrating the GRN with NicheNet priors.
- `ContextResult`: answer to one biological question. It stores ranked ligands, target genes, ligand-target links, and the parameters defining sender, receiver, contrast, and thresholds.

Use accessor functions to extract ordinary R objects:

```r
edges(grn)
tfs(grn)
targets(grn)
regulons(grn)

ligand_target_matrix(network)
ligands(network)
integration_params(network)

ligand_activities(result)
top_ligands(result)
top_targets(result)
ligand_target_links(result)
parameters(result)
as.data.frame(result)
provenance(result)
```

## ccflowR in the Tutorials

ccflowR is used after ligand prioritization to visualize path-level hypotheses. A ccflowR path should be read as a prioritized model-supported route from sender to receiver. It is not direct proof of secretion, receptor binding, spatial contact, or causality. The strongest interpretation comes when the same path is supported by ligand expression, receptor expression, target-gene changes, regulatory potential, abundance information, and external biological evidence.

## Troubleshooting

| Symptom | Likely cause | What to try |
|---------|--------------|-------------|
| Tutorial 01 fails before loading data | Package install or kernel points to the wrong R environment | Run `check_ecosystem(verbose = TRUE)` in the same kernel and verify `ContextNet R` maps to the installed environment |
| Few GRN targets overlap the expression object | Species or gene-symbol mismatch | Harmonize gene naming before `import_grn()`; use official symbols consistently |
| Empty or near-empty ligand-target matrix | GRN is too sparse or all-zero | Inspect `edges(grn)` density and relax upstream GRN inference thresholds |
| Tiny receiver gene set or unstable ranking | DE contrast yields too few genes | Loosen thresholds or pass a curated `geneset_oi` |
| No receptors for top ligands | Receptors are not expressed in the receiver or absent from the prior | Check receptor expression and the NicheNet prior version |
| ccflowR Sankey is empty | Threshold too strict | Lower `threshold_quantile` in `combine_sankey_data()` |
| MultiNicheNet FORK/socket errors | Parallel backend issue | Set `n_cores = 1` |

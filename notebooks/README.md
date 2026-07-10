# ContextNet Tutorials

ContextNet extends NicheNet by replacing the generic gene regulatory network with a context-specific gene regulatory network (GRN). ContextNet is useful when a regulatory network inferred from sc-multiomic data should guide ligand-target predictions instead of a generic prior. The tutorials move from installation verification and object orientation to HNSCC analyses that ask which communication programs may explain epithelial subtype organization and condition-associated signaling programs.

## Tutorial Map

| Tutorial | Question addressed | Data mode | Main output | Suitable use |
|----------|--------------------|-----------|-------------|--------------|
| [01 Getting Started](01-getting-started.ipynb) / [Rmd](01-getting-started.Rmd) | Does ContextNet load correctly, and what do the core inputs and outputs contain? | Included example dataset | `ContextGRN`, `ContextNetwork`, `ContextResult`, expected plots | Installation checks and object-model orientation |
| [02 Standard Workflow](02-standard-workflow.ipynb) / [Rmd](02-standard-workflow.Rmd) | Which communication programs may explain epithelial cell states within one tissue context? | HNSCC Seurat object + context-specific GRN | epithelial-subtype ligand activity, target heatmaps, ccflowR Sankey river paths | Single-tissue ContextNet analysis across epithelial subtypes |
| [03 Multi-Condition Workflow](03-multicondition-workflow.ipynb) / [Rmd](03-multicondition-workflow.Rmd) | Which communication programs change reproducibly across biological samples? | HNSCC replicate-aware MultiNicheNet analysis | prioritization table, abundance information, ccflowR Sankey paths | Sample-aware inference with biological replicates |

RStudio users can open the companion `.Rmd` files in this directory. They contain the same source cells as the notebooks and are meant for interactive execution, not package build-time vignettes.

## Running the Tutorials

- **Installation.** Complete package setup with the repository [installation instructions](../README.md#installation). The tutorials do not install packages; setup cells only load libraries and data.

- **Kernel.** The notebooks use an R kernel named `contextnet_r`, displayed in Jupyter as **ContextNet R**. Register it once with `IRkernel::installspec(name = "contextnet_r", displayname = "ContextNet R")`. If Jupyter also shows an old raw `ir` kernel for the same environment, remove the confusing entry with `jupyter kernelspec remove ir -f` after confirming that `contextnet_r` works.

- **Working directory.** Run a notebook from the repository root or from the `notebooks/` directory. The setup cell normalizes the working directory to the repository root.

- **ccflowR river plots.** [ccflowR](https://github.com/settylab/ccflowR) is optional visualization support for Sankey/river plots. It helps visualize high-priority sender -> ligand -> receptor -> receiver paths after ContextNet has prioritized communication hypotheses. Install it into the same **ContextNet R** environment with `remotes::install_github("settylab/ccflowR", upgrade = "never")`. Without ccflowR, core ContextNet analyses still run, but Sankey/river plot cells should be skipped. See the [ccflowR documentation](https://github.com/settylab/ccflowR#readme) for plotting-specific options.

- **Outputs and figures.** Tutorial 01 is committed with outputs from the included example dataset. Tutorials 02 and 03 are committed with current executed outputs from the HNSCC tutorial analyses so the expected tables and figures are visible. Re-run notebooks whenever inputs, package versions, or plotting code change.

## Data Requirements

Tutorial 01 uses the small example dataset included with the package at `inst/extdata/contextnet_smoke/`. It is synthetic and is intended for installation checks and object interpretation.

Tutorials 02 and 03 use HNSCC tutorial data. Download the data from the [ContextNet tutorial data release](https://github.com/settylab/ContextNet/releases), place the files under `data/`, or set the documented environment variables before running the notebooks.

| Input | File used by the tutorials | Purpose |
|-------|----------------------------|---------|
| Context-specific GRN | `gene_x_tf.csv` | target-gene x transcription-factor regulatory weights from a GRN inference workflow |
| Seurat object | `seurat_obj.rds` | expression matrix plus cell-type, condition, and sample metadata |
| NicheNet L-R prior | `lr_network_human_21122021.rds` | ligand-receptor interactions |
| NicheNet weighted prior | `weighted_networks_nsga2r_final.rds` | signaling and regulatory priors for propagation |
| NicheNet L-T matrix (optional) | `ligand_target_matrix_nsga2r_final.rds` | reference only; not used by ContextNet integration |

The NicheNet priors are distributed with the official [NicheNet v2 tutorials](https://github.com/saeyslab/nichenetr) and can be fetched once with `load_nichenet_data(cache_dir = "data", version = "v2")`. The HNSCC Seurat object and saved MultiNicheNet result are provided as tutorial data rather than bundled inside the R package.

## How to Interpret the Object Model

ContextNet separates the analysis into three object classes:

- `ContextGRN`: regulatory evidence from the profiled biological context. It stores transcription-factor to target-gene edges imported from a GRN inference workflow.
- `ContextNetwork`: reusable propagated signaling model. It stores the context-specific ligand-target matrix created by integrating the GRN with NicheNet priors.
- `ContextResult`: result for one sender/receiver analysis. It stores ranked ligands, target genes, ligand-target links, and the parameters defining sender cells, receiver cells, contrasts, and thresholds.

Accessor functions extract ordinary R objects for inspection, plotting, and export. Tutorial 01 explains these objects with printed output in the [object model summary](01-getting-started.ipynb#Object-Model-Summary).

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

ccflowR is used after ligand prioritization to visualize path-level hypotheses. A ccflowR path should be read as a prioritized model-supported route from sender to receiver. It is not direct proof of secretion, receptor binding, spatial contact, or causality. The strongest interpretation comes when the same path is supported by ligand expression, receptor expression, target-gene changes, regulatory-potential scores, abundance information, and external biological evidence. See the [ccflowR documentation](https://github.com/settylab/ccflowR#readme) for additional plotting examples.

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

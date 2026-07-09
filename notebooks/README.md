# ContextNet Tutorials

The tutorials are organized around biological questions rather than package functions. ContextNet analyses should be read as communication hypotheses: a sender ligand is prioritized when its propagated targets explain a receiver-cell transcriptional program through a context-specific regulatory network.

## Tutorial Map

| Tutorial | Biological question | Statistical mode | Main output | Use when |
|----------|---------------------|------------------|-------------|----------|
| [01 Getting Started](01-getting-started.ipynb) | Which ligands are consistent with a basal epithelial identity program? | Cell type comparison | `ContextResult` with ligand activity and ligand-target links | You are learning the object model and minimal workflow |
| [02 Standard Workflow](02-standard-workflow.ipynb) | Which sender ligands may explain basal, parabasal, suprabasal, and dividing receiver states within contralateral tissue? | Contralateral-only epithelial state comparison | receiver-specific ligand activity, target heatmaps, ccflowR Sankey river paths | You want a single-tissue ContextNet analysis across epithelial receiver states |
| [03 Multi-Condition Workflow](03-multicondition-workflow.ipynb) | Which communication programs change reproducibly across samples? | Replicate-aware MultiNicheNet pseudobulk | prioritization table and ccflowR Sankey paths | You have biological replicates and need sample-aware inference |

## Running the Notebooks

- **Installation.** Install ContextNet first using the repository `README.md` or the micromamba guide. The notebooks do not install packages; they start from biological questions and use setup cells only to load libraries and local data.

- **Kernel.** The notebooks use an R kernel named `contextnet_r`, displayed in Jupyter as **ContextNet R**. Install it once from R with `IRkernel::installspec(name = "contextnet_r", displayname = "ContextNet R")` (see the [micromamba installation guide](../vignettes/installation-micromamba.md)). If your kernel has a different name, either select it in Jupyter or override it on the command line:

  ```bash
  jupyter nbconvert --to notebook --execute \
      --ExecutePreprocessor.kernel_name=contextnet_r \
      notebooks/01-getting-started.ipynb
  ```

- **Working directory.** Run a notebook from the repository root or from the `notebooks/` directory. The first setup cell normalizes the working directory to the repository root, so all data paths are relative and portable.

- **ccflowR river plots.** [ccflowR](https://github.com/settylab/ccflowR) is required for the committed Sankey river plots in Tutorials 02 and 03. Install it into the same **ContextNet R** kernel environment before running the notebooks; the micromamba guide shows the private-repo install step.

- **Outputs and figures.** The notebooks are committed with current executed outputs from the local HNSCC example data so readers can see the expected tables and figures. Re-run a notebook locally when inputs, package versions, or plotting code change.

## Data Requirements and Example-Data Strategy

These are **full local-data tutorials**: they read inputs that are present locally rather than shipping a bundled toy dataset. To run them as written you need, under `data/`:

| Input | File (as used in the tutorials) | Purpose |
|-------|----------------------------------|---------|
| Context-specific GRN | `gene_x_tf.csv` | gene × TF regulatory weights from a GRN inference workflow |
| Seurat object | `seurat_obj.rds` | expression + cell-type / condition / sample metadata |
| NicheNet L-R prior | `lr_network_human_21122021.rds` | ligand-receptor interactions |
| NicheNet weighted prior | `weighted_networks_nsga2r_final.rds` | signaling + regulatory priors for propagation |
| NicheNet L-T matrix (optional) | `ligand_target_matrix_nsga2r_final.rds` | reference only; not used by integration |

The NicheNet priors are distributed with the official NicheNet v2 tutorials ([Zenodo 7074291](https://zenodo.org/record/7074291)) and can be fetched once with `load_nichenet_data(cache_dir = "data", version = "v2")`. The GRN and Seurat object are your own data; the tutorials use a head-and-neck squamous cell carcinoma (HNSCC) dataset as the worked example.

**Why no bundled example data?** A realistic run requires a multi-cell-type Seurat object plus the full NicheNet priors — the local `data/` directory these tutorials read is **~1.8 GB in total** (the Seurat object alone is ~1.2 GB; the optional NicheNet ligand-target matrix is >250 MB), far too large to ship in the package. Provision this space before you start. If you want a fast smoke test, subset your own Seurat object before running, for example:

```r
# Small, fast example: keep a few cell types and down-sample cells.
keep <- WhichCells(seurat_obj, expression = Fine_Cell_Types %in%
                   c("basal", "parabasal", "suprabasal", "Fibroblast", "Myeloid"))
small <- subset(seurat_obj, cells = sample(keep, min(8000, length(keep))))
saveRDS(small, "data/seurat_obj_small.rds")
```

Point `SEURAT_FILE` at the subset. Ligand activity rankings from a small subset are for pipeline testing only, not biological interpretation.

## How to Interpret the Object Model

ContextNet separates the analysis into three object classes:

- `ContextGRN`: the regulatory evidence from the profiled biological context. It stores TF-target edges imported from a GRN inference workflow.
- `ContextNetwork`: the reusable propagated signaling model. It stores the context-specific ligand-target matrix created by integrating the GRN with NicheNet priors.
- `ContextResult`: the answer to one specific biological question. It stores ranked ligands, target genes, ligand-target links, and the parameters that define the sender, receiver, contrast, and thresholds.

Use accessor functions to extract data:

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

The extracted values are ordinary R objects and are the recommended inputs for custom plots, enrichment analysis, external validation, and direct use of `nichenetr` or `multinichenetr`.

## ccflowR in the Tutorials

ccflowR is used after ligand prioritization to visualize path-level hypotheses. Tutorial 02 maps receiver-specific `ContextResult` objects into Sankey rows:

```r
sankey_raw <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)
ccflowR::sankey_flow(sankey_data, cell_colors)
```

Tutorial 03 follows the supplied MultiNicheNet reference analysis and plots the top rows of `group_prioritization_tbl` directly with `ccflowR::sankey_plot()`.

Read a ccflowR path as a prioritized model-supported route from sender to receiver. It is not direct proof of secretion, receptor binding, or spatial contact. The strongest interpretation comes when the same path is supported by ligand expression, receptor expression, target-gene changes, regulatory potential, and external biological evidence.

## Troubleshooting

| Symptom | Likely cause | What to try |
|---------|--------------|-------------|
| Few GRN targets overlap the expression object | Species or gene-symbol mismatch (e.g. mouse vs human, Ensembl IDs vs symbols) | Harmonize gene naming before `import_grn()`; the notebooks warn when overlap is low |
| Empty or near-empty ligand-target matrix | GRN is too sparse or all-zero | Inspect `edges(grn)` density; relax GRN inference thresholds upstream |
| Tiny receiver gene set / unstable ranking | DE contrast yields few genes | Loosen `logfc_threshold` / `pval_threshold`, or pass a curated `geneset_oi` |
| No receptors for top ligands | Receptors not expressed in the receiver, or absent from the L-R prior | Check receptor expression; confirm the correct NicheNet prior version |
| ccflowR Sankey is empty | Threshold too strict | Lower `threshold_quantile` in `combine_sankey_data()` |
| MultiNicheNet FORK/socket errors | Parallel backend issues on the platform | Set `n_cores = 1` (Tutorial 03 default) |

## Remaining Documentation Gaps to Address Next

- Add a dedicated GRN-input tutorial covering contextnet-pipeline, SCENIC+, Pando, and custom edge tables.
- Add a comparison tutorial that runs the same contrast with the default NicheNet matrix and the ContextNet matrix to show what biological specificity changes.
- Consider a packaged minimal example object (`inst/extdata/`) once a sufficiently small, shareable dataset is available.

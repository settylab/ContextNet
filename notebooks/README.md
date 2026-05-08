# ContextNet Tutorials

The tutorials are organized around biological questions rather than package functions. ContextNet analyses should be read as communication hypotheses: a sender ligand is prioritized when its propagated targets explain a receiver-cell transcriptional program through a context-specific regulatory network.

## Tutorial Map

| Tutorial | Biological question | Statistical mode | Main output | Use when |
|----------|---------------------|------------------|-------------|----------|
| [01 Getting Started](01-getting-started.ipynb) | Which ligands are consistent with a basal epithelial identity program? | Cell type comparison | `ContextResult` with ligand activity and ligand-target links | You are learning the object model and minimal workflow |
| [02 Standard Workflow](02-standard-workflow.ipynb) | Which sender ligands may explain basal-cell changes between C and N tissue? | Within-receiver condition comparison | ligand activity, target heatmaps, ccflowR Sankey paths | You have a single contrast and want full interpretation/visualization |
| [03 Multi-Condition Workflow](03-multicondition-workflow.ipynb) | Which communication programs change reproducibly across samples? | Replicate-aware MultiNicheNet pseudobulk | prioritization table and ccflowR Sankey paths | You have biological replicates and need sample-aware inference |

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

ccflowR is used after ligand prioritization to visualize path-level hypotheses:

```r
sankey_raw <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)
ccflowR::sankey_flow(sankey_data, cell_colors)
```

Read a ccflowR path as a prioritized model-supported route from sender to receiver. It is not direct proof of secretion, receptor binding, or spatial contact. The strongest interpretation comes when the same path is supported by ligand expression, receptor expression, target-gene changes, regulatory potential, and external biological evidence.

## Remaining Documentation Gaps to Address Next

- Add small bundled example data so the notebooks can run without private local paths.
- Add a dedicated GRN-input tutorial covering contextnet-pipeline, SCENIC+, Pando, and custom edge tables.
- Add a troubleshooting section for gene-name mismatches, sparse GRNs, missing receptors, and empty ccflowR outputs.
- Add a comparison tutorial that runs the same contrast with the default NicheNet matrix and the ContextNet matrix to show what biological specificity changes.

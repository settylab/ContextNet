# ContextNet

**scMultiomics-Enhanced Cell-Cell Communication Inference**

ContextNet extends [NicheNet](https://github.com/saeyslab/nichenetr) by replacing its generic gene regulatory network with a **context-specific GRN** inferred from single-cell multiomics data (scRNA-seq + scATAC-seq). This grounds ligand-target predictions in experimentally observed regulatory relationships specific to your biological system.

## Overview

NicheNet models intercellular communication by linking ligands to target genes through three network layers: ligand-receptor interactions, intracellular signaling, and gene regulation. Its gene regulatory network derives from aggregated prior knowledge across diverse biological contexts. ContextNet replaces this generic GRN with one inferred directly from your data, while preserving NicheNet's curated signaling and ligand-receptor networks.

The framework has three components:

| Component | Description |
|-----------|-------------|
| **GRN inference workflow** | Snakemake pipeline for constructing context-specific TF-target relationships from matched scRNA-seq + scATAC-seq |
| **ContextNet** (this package) | R package for GRN integration with NicheNet and ligand activity prediction |
| **[ccflowR](https://github.com/settylab/ccflowR)** | R package for Sankey flow visualization of cell-cell communication |

## Installation

```r
# Install ContextNet (includes NicheNet and MultiNicheNet)
remotes::install_github("settylab/ContextNet", dependencies = TRUE)
```

For a fully reproducible environment using **micromamba** (recommended for development), see the [Micromamba Installation Guide](vignettes/installation-micromamba.md).

### Verify installation

```r
library(ContextNet)
check_ecosystem(verbose = TRUE)
```

## Quick Start

```r
library(ContextNet)

# 1. Import your context-specific GRN (gene x TF matrix, CSV)
grn <- import_grn("gene_x_tf.csv")

# 2. Integrate with NicheNet networks
network <- integrate_networks(grn, nichenet_database = "data/")

# 3. Predict ligand activity
result <- predict_ligand_activity_nichenet(
    object = seurat_obj,
    network = network,
    sender = "Myeloid",
    receiver = "Fibroblast",
    group_col = "cell_type",
    ident_col = "condition",
    ident_1 = "Disease",
    ident_2 = "Control",
    nichenet_database = "data/"
)

# 4. Explore results
top_ligands(result)
ligand_activities(result)
```

### Sankey visualization with ccflowR

```r
# Prepare per-receiver data
data_recv <- prepare_for_ccflowr(result, seurat_obj, network)

# Combine, scale globally, threshold
sankey_data <- combine_sankey_data(data_recv, threshold_quantile = 0.95)

# Visualize (requires ccflowR)
library(ccflowR)
ccflowR::sankey_flow(sankey_data, cell_colors)
```

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [01 - Getting Started](notebooks/01-getting-started.ipynb) | Concepts, installation, and quick 4-step workflow |
| [02 - Standard Workflow](notebooks/02-standard-workflow.ipynb) | Full NicheNet workflow with visualization and export |
| [03 - Multi-Condition Workflow](notebooks/03-multicondition-workflow.ipynb) | Differential analysis with MultiNicheNet and biological replicates |

## Key Functions

### GRN Import

```r
grn <- import_grn("gene_x_tf.csv")  # Gene x TF matrix (CSV)
```

### Network Integration

```r
network <- integrate_networks(grn, nichenet_database = "data/", mode = "replace")
```

Three modes: `replace` (default, most context-specific), `augment`, `weighted`.

### Communication Analysis

```r
# Standard NicheNet (single condition or no replicates)
result <- predict_ligand_activity_nichenet(
    object, network, sender, receiver,
    ident_col = "condition", ident_1 = "Disease", ident_2 = "Control"
)

# MultiNicheNet (multi-condition with replicates)
result <- predict_ligand_activity_multinichenet(
    object, network, sender, receiver,
    ident_col = "condition", ident_1 = "Disease", ident_2 = "Control",
    sample_col = "sample"  # Required
)
```

### Visualization Export

```r
# Sankey flow data for ccflowR
data <- prepare_for_ccflowr(result, seurat_obj, network)
sankey <- combine_sankey_data(data_recv1, data_recv2, threshold_quantile = 0.95)

# Cytoscape
export_to_cytoscape(network, "network.graphml")

# Save/load ContextNet objects
save_contextnet(result, "analysis.cnet")
result <- load_contextnet("analysis.cnet")
```

## Related Packages

- **[nichenetr](https://github.com/saeyslab/nichenetr)** — NicheNet: ligand-target modeling via prior knowledge networks
- **[multinichenetr](https://github.com/saeyslab/multinichenetr)** — MultiNicheNet: multi-condition differential communication analysis
- **[ccflowR](https://github.com/settylab/ccflowR)** — Sankey flow visualization for cell-cell communication
- **[Seurat](https://satijalab.org/seurat/)** — Single-cell RNA-seq analysis toolkit

## Citation

If you use ContextNet in your research, please cite:

> [Citation to be added upon publication]

ContextNet builds on NicheNet and MultiNicheNet:

> Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. *Nature Methods* 17, 159-162 (2020).

> Browaeys R et al. MultiNicheNet: a flexible framework for cell-cell communication analysis from multi-sample multi-condition single-cell transcriptomics data. *Nature Methods* (2024).

## License

MIT

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

The recommended path is the bundled micromamba environment, which produces a reproducible install of R, every Bioconductor and CRAN dep, and the GitHub-only peers (`nichenetr`, `multinichenetr`, `ccflowR`):

```bash
# 1. Create the conda environment (R + system libs + most CRAN/Bioconductor deps)
micromamba create -n contextnet -f environment.yml -y

# 2. Bootstrap the GitHub-only deps and ContextNet itself (idempotent)
#    ccflowR is in a private repository; export a PAT with `repo` scope.
export GITHUB_PAT=$(gh auth token)         # or paste a classic PAT
micromamba run -n contextnet Rscript inst/bootstrap.R
```

The bootstrap script installs the GitHub peers at audited pins, including `settylab/ccflowR` at commit `6056c0fe1d1cecfc68aa059bacfbe6bd16725231`. Without `GITHUB_PAT` it skips `ccflowR` (the Sankey visualization in tutorials 02/03 will be unavailable) and continues. See the [Micromamba Installation Guide](vignettes/installation-micromamba.md) for the full reference, manual install steps, and troubleshooting.

### Lightweight install (CRAN-only path)

If you don't need the multinichenetr Bioconductor stack and have a working R toolchain, you can install ContextNet alone with `remotes`. You will still need to install `ccflowR` manually if you want Sankey visualization:

```r
# In an R session with `remotes` and a working compiler toolchain
remotes::install_github("settylab/ContextNet", dependencies = TRUE)

# Optional: ccflowR for Sankey visualization (private repo — set GITHUB_PAT first)
Sys.setenv(GITHUB_PAT = "<your-token>")
remotes::install_github("settylab/ccflowR",
                        ref = "6056c0fe1d1cecfc68aa059bacfbe6bd16725231")
```

### Verify installation

```r
library(ContextNet)
check_ecosystem(verbose = TRUE)
```

### Known installation issues

| Symptom | Cause / Fix |
|---|---|
| `remotes::install_github('settylab/ccflowR')` returns 404 | Private repository. Set `GITHUB_PAT` (`$(gh auth token)` works if you have `gh` logged in) and retry. ContextNet works without ccflowR; you only lose the tutorial Sankey visualization. |
| `nloptr` fails to compile from source | The micromamba env ships `nlopt` + `r-nloptr` pre-built. If you're outside the env: `micromamba install -n contextnet -c conda-forge nlopt r-nloptr r-lme4 -y`. |
| `hdf5r` configure: `hdf5 does not seem to be installed` | The micromamba env ships `r-hdf5r` pre-built with bundled HDF5. Recreate the env from `environment.yml`. Outside the env: install your distro's `libhdf5-dev` / `hdf5-devel`, or `brew install hdf5`. |
| `failed to lock directory` after a killed install | Stale `00LOCK-<pkg>` directories. Re-run `Rscript inst/bootstrap.R` (it clears them automatically) or `rm -rf $CONDA_PREFIX/lib/R/library/00LOCK-*`. |
| `repo.anaconda.com` unreachable | The shipped `environment.yml` only uses `conda-forge` and `bioconda`. If a global `~/.condarc` adds the default channel, pass `--override-channels` to `micromamba create`. |
| macOS: "Operation not permitted" on the working tree | Extended attributes blocking the build: `xattr -cr /path/to/ContextNet`. |

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

ContextNet pairs with [ccflowR](https://github.com/settylab/ccflowR) to render cell-cell communication as a Sankey diagram (sender → ligand → receptor → target → receiver).

```r
library(ccflowR)

# 1. Build per-receiver Sankey data: composite score from ligand expression and predicted activity
sankey_raw <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)

# 2. Combine across receivers, scale globally, threshold to top quantile
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)

# 3. Render
cell_colors <- c(Myeloid = "#E41A1C", Fibroblast = "#377EB8", Epithelial = "#4DAF4A")
sankey_flow(sankey_data, cell_colors)
```

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [01 - Getting Started](notebooks/01-getting-started.ipynb) | Concepts, installation, and quick 4-step workflow |
| [02 - Standard Workflow](notebooks/02-standard-workflow.ipynb) | Full NicheNet workflow with Sankey visualization |
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

### Visualization

```r
# Build Sankey data and render with ccflowR
sankey_raw  <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)
ccflowR::sankey_flow(sankey_data, cell_colors)
```

### Save / Load

```r
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

GPL (>= 3). See [LICENSE.md](LICENSE.md) for the full license text.

ContextNet imports [nichenetr](https://github.com/saeyslab/nichenetr) (GPL-3, Saeys lab, Ghent University), so the package is distributed under a GPL-compatible license.

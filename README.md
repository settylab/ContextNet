# ContextNet

**scMultiomics-Enhanced Cell-Cell Communication Inference**

ContextNet extends [NicheNet](https://github.com/saeyslab/nichenetr) by replacing its generic gene regulatory network with a **context-specific gene regulatory network (GRN)** inferred from single-cell multiomics data (scRNA-seq + scATAC-seq). This grounds ligand-target predictions in experimentally observed regulatory relationships specific to your biological system.

## Overview

NicheNet models intercellular communication by linking ligands to target genes through three network layers: ligand-receptor interactions, intracellular signaling, and gene regulation. Its gene regulatory network derives from aggregated prior knowledge across diverse biological contexts. ContextNet replaces this generic GRN with one inferred directly from your data, while preserving NicheNet's curated signaling and ligand-receptor networks.

The framework has three components:

| Component | Description |
|-----------|-------------|
| **GRN inference workflow** | Snakemake pipeline for constructing context-specific transcription-factor to target-gene relationships from matched scRNA-seq + scATAC-seq. Public pipeline link: <https://github.com/settylab/contextnet-pipeline> |
| **ContextNet** (this package) | R package for GRN integration with NicheNet and ligand activity prediction |
| **[ccflowR](https://github.com/settylab/ccflowR)** | Optional R package for Sankey/river visualization of sender-ligand-receptor-receiver communication paths |

## Installation

The recommended path is the bundled micromamba environment. It installs R, compiled system libraries, CRAN packages, Bioconductor packages, Jupyter, and IRkernel from conda-forge/bioconda where possible, then uses `inst/bootstrap.R` for GitHub-only R packages.

```bash
# 1. Create the conda environment
micromamba create -n contextnet -f environment.yml -y

# 2. Install GitHub-only dependencies and ContextNet itself
micromamba run -n contextnet Rscript inst/bootstrap.R

# 3. Register the notebook kernel with an intuitive name
micromamba run -n contextnet R -e "IRkernel::installspec(name = 'contextnet_r', displayname = 'ContextNet R')"
```

The bootstrap installs required GitHub packages such as `nichenetr` and `multinichenetr`, then installs ContextNet from the local checkout. It does **not** require private GitHub access. This is important for publication use: a reader outside the lab should be able to install the core package without a settylab token.

`ccflowR` is optional visualization support for Sankey/river plots. It turns ContextNet or MultiNicheNet result tables into path-level diagrams that connect sender cell type, ligand, receptor, and receiver cell type. This is useful when ranked ligand tables are too dense and the reader needs to inspect the strongest cell-cell communication routes. Once `ccflowR` is public, install it into the same R environment with:

```bash
micromamba run -n contextnet R -e "remotes::install_github('settylab/ccflowR', upgrade = 'never')"
```

or opt in during bootstrap:

```bash
INSTALL_CCFLOWR=true micromamba run -n contextnet Rscript inst/bootstrap.R
```

If `ccflowR` is not installed, ContextNet still installs and the core analyses run; only Sankey/river plot cells are unavailable. See the [Micromamba Installation Guide](vignettes/installation-micromamba.md) for manual install steps, kernel registration, and troubleshooting.

### Lightweight install

If you already have a working R toolchain and only need ContextNet core functionality, install with hard dependencies first. Avoid `dependencies = TRUE` unless optional packages such as `ccflowR`, `Seurat`, and the full MultiNicheNet stack are already available.

```r
# In an R session with remotes installed
remotes::install_github("settylab/ContextNet", dependencies = NA)

# Optional Sankey/river-plot support once ccflowR is public
remotes::install_github("settylab/ccflowR", upgrade = "never")
```

### Verify installation

```r
library(ContextNet)
check_ecosystem(verbose = TRUE)
```

For a fast functional test, run Tutorial 01. It uses the packaged synthetic smoke-test fixture under `inst/extdata/contextnet_smoke/`, so it does not require lab-private Seurat or AnnData files.

### Known installation issues

| Symptom | Cause / Fix |
|---|---|
| `ccflowR` is not available | ContextNet core analysis does not require it. Install it from `settylab/ccflowR` once the repository is public; until then, skip Sankey/river plots. |
| Jupyter shows a raw `ir` kernel name | Register the intended kernel with `IRkernel::installspec(name = 'contextnet_r', displayname = 'ContextNet R')`. If an old `ir` kernelspec points to the same env and causes confusion, remove it with `jupyter kernelspec remove ir -f`. |
| `nloptr` fails to compile from source | Use the micromamba env, which ships `nlopt` and `r-nloptr` pre-built. Outside the env, install system `nlopt` first. |
| `hdf5r` configure: `hdf5 does not seem to be installed` | Use the micromamba env, which ships `r-hdf5r` with HDF5. Outside the env, install your platform HDF5 development libraries. |
| `failed to lock directory` after a killed install | Remove stale `00LOCK-*` directories or re-run `Rscript inst/bootstrap.R`, which clears them. |
| `repo.anaconda.com` unreachable | The shipped `environment.yml` uses only `conda-forge` and `bioconda`. If global conda config adds defaults, pass `--override-channels`. |
| macOS: "Operation not permitted" on the working tree | Clear extended attributes with `xattr -cr /path/to/ContextNet`. |

## Quick Start

```r
library(ContextNet)

# 1. Import your context-specific gene regulatory network (GRN)
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

ContextNet pairs with [ccflowR](https://github.com/settylab/ccflowR) to render cell-cell communication as a Sankey river diagram: sender cell type -> ligand -> receptor -> receiver cell type. ccflowR is helpful when a ranked ligand table is too dense and you need to inspect the strongest path-level communication hypotheses. `prepare_for_ccflowr()` and `combine_sankey_data()` build the Sankey data table from base ContextNet, and `ccflowR::sankey_plot()` renders the river plot used in Tutorials 02 and 03.

```r
library(ccflowR)

# 1. Build per-receiver Sankey data: composite score from ligand expression and predicted activity
sankey_raw <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)

# 2. Combine across receivers, scale globally, threshold to top quantile
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)

# 3. Render
cell_colors <- c(Myeloid = "#E41A1C", Fibroblast = "#377EB8", Epithelial = "#4DAF4A")
ccflowR::sankey_plot(sankey_data, colors = cell_colors)
```

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [01 - Getting Started](notebooks/01-getting-started.ipynb) ([Rmd](notebooks/01-getting-started.Rmd)) | Public smoke test for installation, input files, output objects, and expected figures |
| [02 - Standard Workflow](notebooks/02-standard-workflow.ipynb) ([Rmd](notebooks/02-standard-workflow.Rmd)) | Contralateral-only receiver-state workflow with ccflowR Sankey river visualization |
| [03 - Multi-Condition Workflow](notebooks/03-multicondition-workflow.ipynb) ([Rmd](notebooks/03-multicondition-workflow.Rmd)) | Differential analysis with MultiNicheNet and biological replicates |

Installation is a separate setup step covered above and in the [Micromamba Installation Guide](vignettes/installation-micromamba.md). The notebooks assume the package is already installed, then walk through biological questions using local data. Tutorial 01 uses a packaged synthetic fixture. Tutorials 02 and 03 read a context-specific gene regulatory network (GRN), a Seurat object, and the NicheNet priors from local data paths and do not download anything while running. See [notebooks/README.md](notebooks/README.md) for data requirements, kernel setup, expected outputs, and the R Markdown companion files.

## Key Functions

### Gene Regulatory Network (GRN) Import

```r
# Target-gene x transcription-factor matrix from the GRN inference pipeline
grn_file <- "gene_x_tf.csv"
grn <- import_grn(grn_file)
```

A context-specific gene regulatory network (GRN) is the core input that makes ContextNet different from default NicheNet. It stores transcription-factor to target-gene regulatory weights inferred from the profiled biological context. Use the ContextNet GRN inference pipeline (<https://github.com/settylab/contextnet-pipeline>) or provide an equivalent matrix or edge table from another regulatory-inference workflow.

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
result <- run_multinichenet(
    object, network,
    sender_receiver_pairs = sender_receiver_pairs,
    condition_col = "condition",
    contrast = c(condition_oi = "Disease", condition_ref = "Control"),
    sample_col = "sample"  # Required
)
```

### Visualization

```r
# Build Sankey data and render with ccflowR for a ContextResult
sankey_raw  <- prepare_for_ccflowr(result, expression_data = seurat_obj, network = network)
sankey_data <- combine_sankey_data(sankey_raw, threshold_quantile = 0.95)
ccflowR::sankey_plot(sankey_data, colors = cell_colors)
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

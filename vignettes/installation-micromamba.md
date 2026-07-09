# Installation with Micromamba

This guide covers how to install the ContextNet R environment with micromamba. It is intentionally limited to installation, kernel registration, verification, and debugging. Biological questions and interpretation start in the tutorial notebooks and R Markdown files after the environment works.

## Why Micromamba?

ContextNet depends on R, Bioconductor packages, GitHub-only R packages, and compiled system libraries such as `nlopt`, `Cairo`, and HDF5. Micromamba installs most of these as pre-built binaries from conda-forge and bioconda, which avoids common source-build failures on Linux and macOS.

## Prerequisites

Install micromamba, then restart your shell.

```bash
# macOS with Homebrew
brew install micromamba

# Linux or macOS with the installer
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Initialize your shell, then restart the terminal
micromamba shell init -s bash  # or zsh/fish

micromamba --version
```

Windows users should follow the official micromamba installation documentation.

## Create the Environment

Run from the ContextNet repository root:

```bash
micromamba create -n contextnet -f environment.yml -y
```

This installs:

- R and the R build toolchain,
- required CRAN packages,
- required Bioconductor packages available from conda,
- system libraries needed by R packages,
- Jupyter and IRkernel for notebooks,
- `rmarkdown` and `knitr` for RStudio/R Markdown users.

If a global conda configuration adds the Anaconda default channel and your network blocks it, override channels explicitly:

```bash
micromamba create -n contextnet -f environment.yml --override-channels -y
```

## Install ContextNet and Required GitHub Packages

After creating the environment, run the bootstrap script. It installs required GitHub packages (`nichenetr`, `multinichenetr`), Bioconductor dependencies needed by MultiNicheNet, and ContextNet itself from the local checkout.

```bash
micromamba run -n contextnet Rscript inst/bootstrap.R
```

The bootstrap is idempotent. It skips packages that already satisfy the required version, clears stale `00LOCK-*` directories from interrupted installs, and then runs `check_ecosystem()` so you can see which optional components are available.

The bootstrap uses hard dependencies for the ContextNet install. That keeps installation reliable for external users because optional visualization packages do not block the core package install.

## Optional: Install ccflowR for Sankey/River Plots

`ccflowR` is optional visualization support. It renders sender -> ligand -> receptor -> receiver communication paths as Sankey/river plots, which is helpful when the result table contains many plausible cell-cell communication routes and you need to inspect the highest-priority paths. ContextNet core analysis does not require `ccflowR`; only Sankey/river plot cells require it.

Once the `settylab/ccflowR` repository is public, install it into the same environment:

```bash
micromamba run -n contextnet R -e "remotes::install_github('settylab/ccflowR', upgrade = 'never')"
```

Or opt in during bootstrap:

```bash
INSTALL_CCFLOWR=true micromamba run -n contextnet Rscript inst/bootstrap.R
```

If `ccflowR` is unavailable, do not use a lab-private token as part of the public installation path. Install ContextNet core, run Tutorial 01, and skip Sankey/river plots until `ccflowR` is public.

## Register the Jupyter Kernel

Register a clearly named R kernel so users do not see a generic `ir` entry and wonder which environment is active:

```bash
micromamba run -n contextnet R -e "IRkernel::installspec(name = 'contextnet_r', displayname = 'ContextNet R')"
```

The notebooks are committed with this kernelspec:

```text
name: contextnet_r
display_name: ContextNet R
language: R
```

If Jupyter still shows an old raw `ir` kernel for the same environment, remove it after confirming `ContextNet R` works:

```bash
micromamba run -n contextnet jupyter kernelspec list
micromamba run -n contextnet jupyter kernelspec remove ir -f
```

Run a notebook from the command line with:

```bash
micromamba run -n contextnet jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=contextnet_r \
    notebooks/01-getting-started.ipynb
```

## RStudio and R Markdown

RStudio users can open the companion `.Rmd` files in `notebooks/`. Start RStudio from the activated environment or configure RStudio to use the micromamba R binary:

```bash
micromamba activate contextnet
which R
```

Then run chunks interactively. The `.Rmd` files are tutorial companions, not package-build vignettes, because Tutorials 02 and 03 depend on full local analysis data that is too large to ship with the R package.

## Verify the Installation

First verify the package ecosystem:

```bash
micromamba run -n contextnet R -e "library(ContextNet); check_ecosystem(verbose = TRUE)"
```

Then run Tutorial 01. It uses the packaged synthetic smoke-test fixture under `inst/extdata/contextnet_smoke/` and should not require any lab-private Seurat, AnnData, or NicheNet files. A successful run should produce:

- a `ContextGRN` with 4 regulatory edges,
- a `ContextNetwork` with a small ligand-target matrix,
- a `ContextResult` where `TGFB1` ranks above `FGF2`,
- static ligand activity and ligand-target support plots,
- a Sankey/river plot only if `ccflowR` is installed.

## Manual Install Steps

Use these commands if you need to debug the bootstrap one step at a time.

```bash
# Activate the environment
micromamba activate contextnet

# Required installers
R -e "install.packages('remotes', repos = 'https://cloud.r-project.org')"
R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"

# Bioconductor packages needed by MultiNicheNet
R -e "BiocManager::install(c('muscat', 'scater', 'sva', 'scran'), ask = FALSE, update = FALSE)"

# Required GitHub packages
R -e "remotes::install_github('saeyslab/nichenetr', upgrade = 'never')"
R -e "remotes::install_github('saeyslab/multinichenetr', upgrade = 'never')"

# Optional performance helper for faster Seurat Wilcoxon tests
R -e "remotes::install_github('immunogenomics/presto', upgrade = 'never')"

# ContextNet from the local checkout
R -e "remotes::install_local('.', dependencies = NA, upgrade = 'never')"
```

## Common Issues

| Problem | Solution |
|---------|----------|
| `BiocManager` not found | Recreate the env from `environment.yml` or run `micromamba install -n contextnet -c conda-forge r-biocmanager -y`. |
| `hdf5r` configure error: "hdf5 does not seem to be installed" | Use the micromamba env, which includes `r-hdf5r`; outside the env, install HDF5 development libraries first. |
| `nloptr` fails to compile | Use the micromamba env, which includes `nlopt` and `r-nloptr`; outside the env, install system `nlopt` first. |
| `ccflowR` cannot be installed | ContextNet core still works. Install `ccflowR` from the public GitHub repo once available, then rerun only the Sankey/river plot cells. |
| `failed to lock directory` after an interrupted install | Re-run `Rscript inst/bootstrap.R`, or manually remove `00LOCK-*` directories from the first `.libPaths()` entry. |
| `repo.anaconda.com` unreachable | Use `--override-channels`; this repo's `environment.yml` only requires conda-forge and bioconda. |
| Jupyter uses the wrong R | Register and select `ContextNet R`; remove the confusing `ir` kernelspec if it points to the same environment. |
| RStudio uses a different R | Start RStudio from the activated environment or configure RStudio to use `micromamba run -n contextnet which R`. |

## Quick Reference

```bash
micromamba create -n contextnet -f environment.yml -y
micromamba run -n contextnet Rscript inst/bootstrap.R
micromamba run -n contextnet R -e "IRkernel::installspec(name = 'contextnet_r', displayname = 'ContextNet R')"
micromamba run -n contextnet R -e "library(ContextNet); check_ecosystem(verbose = TRUE)"
micromamba run -n contextnet jupyter nbconvert --to notebook --execute --ExecutePreprocessor.kernel_name=contextnet_r notebooks/01-getting-started.ipynb
```

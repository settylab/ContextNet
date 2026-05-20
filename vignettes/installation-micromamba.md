# Installation with Micromamba

This guide covers how to install and manage the ContextNet development environment using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html), a fast, lightweight package manager compatible with conda environments.

## Why Micromamba?

ContextNet depends on R, Bioconductor packages, and compiled libraries (e.g., `nlopt`, `Cairo`) that are difficult to install from source. Micromamba resolves these as pre-built binaries from conda-forge and bioconda, avoiding common compilation failures on macOS and Linux.

## Prerequisites

### Install micromamba

**macOS (Homebrew):**

```bash
brew install micromamba
```

**Linux / macOS (curl):**

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

**Windows:** See [official docs](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

After installation, initialize your shell:

```bash
micromamba shell init -s bash  # or zsh, fish
# Restart your terminal
```

Verify:

```bash
micromamba --version
```

## Create the Environment

ContextNet provides an `environment.yml` that pins all required R and system dependencies.

```bash
# From the ContextNet project root
micromamba create -n contextnet -f environment.yml -y
```

This installs:
- R (>= 4.1.0) with all required and suggested CRAN packages
- Bioconductor packages (limma, edgeR, locfit)
- System libraries (nlopt, Cairo, gdtools)
- Jupyter + IRkernel for running tutorial notebooks

### Troubleshooting: Channel errors

If you see errors resolving `repo.anaconda.com`, your network may block the Anaconda default channel. The provided `environment.yml` already restricts channels to `conda-forge` and `bioconda`. If you have a global `.condarc` or `~/.mambarc` that adds default channels, override it:

```bash
micromamba create -n contextnet -f environment.yml --override-channels -y
```

## Install Additional Dependencies

Some packages are not available via conda and must be installed from R or GitHub. ContextNet ships a bootstrap script that handles this in one idempotent step.

### Recommended: one-shot bootstrap

The package ships an idempotent bootstrap script that installs the Bioconductor
and GitHub-only dependencies that `environment.yml` cannot supply, then installs
ContextNet itself and runs `check_ecosystem()`. Run it once after creating the
environment:

```bash
micromamba activate contextnet
Rscript inst/bootstrap.R
```

It is safe to re-run after upgrading any single dependency; packages already at a
compatible version are skipped. If you prefer to install the pieces by hand, the
manual steps below do the same thing.

### Activate the environment

```bash
micromamba activate contextnet
```

Or, to run commands without activating:

```bash
micromamba run -n contextnet R -e "<command>"
```

### Run the bootstrap script

From the ContextNet project root:

```bash
# Optional: ccflowR (Sankey visualization) is in a private repository.
# Export a GitHub PAT with `repo` scope to install it. Without a PAT,
# the bootstrap skips ccflowR and continues — ContextNet itself works,
# only the Sankey visualization in tutorials 02/03 will be unavailable.
export GITHUB_PAT=$(gh auth token)   # or paste a classic PAT

micromamba run -n contextnet Rscript inst/bootstrap.R
```

The bootstrap script:

1. Ensures `remotes` and `BiocManager` are installed.
2. Installs the Bioconductor deps `multinichenetr` requires that bioconda does not ship in compatible versions: `muscat`, `scater`, `sva`, `scran`.
3. Installs the GitHub-only packages at pinned versions:
   - `saeyslab/nichenetr` (>= 2.2.1)
   - `saeyslab/multinichenetr` (>= 2.1.0)
   - `settylab/ccflowR` at the audited commit SHA (skipped if `GITHUB_PAT` is unset)
4. Installs `ContextNet` itself from the local source.
5. Runs `check_ecosystem()` to verify every component loads.

The script is idempotent — re-run it after upgrading a single dependency and it will reinstall only what changed. It also clears stale `00LOCK-*` directories from a previous interrupted install, so killing it mid-flight is recoverable.

### Manual install (advanced / offline / debugging)

If you need to inspect each step or you cannot run the bootstrap, the equivalent commands are:

```bash
# Bioconductor deps (multinichenetr's compatible-version requirements)
micromamba run -n contextnet R -e "
  BiocManager::install(
    c('muscat', 'scater', 'sva', 'scran'),
    ask = FALSE, update = FALSE
  )
"

# GitHub packages — ccflowR at audited commit SHA, requires GITHUB_PAT
micromamba run -n contextnet R -e "
  remotes::install_github('saeyslab/nichenetr', upgrade = 'never')
  remotes::install_github('saeyslab/multinichenetr', upgrade = 'never')
  remotes::install_github('settylab/ccflowR',
                          ref = '6056c0fe1d1cecfc68aa059bacfbe6bd16725231',
                          upgrade = 'never')
"

# ContextNet itself (from a local checkout)
micromamba run -n contextnet R -e "
  remotes::install_local('.', dependencies = TRUE, upgrade = 'never')
"
```

### Register the Jupyter R kernel (for the tutorials)

The tutorial notebooks reference the standard IRkernel under the name `ir`.
Register it once so Jupyter (and `nbconvert`) can find it:

```bash
micromamba run -n contextnet R -e "IRkernel::installspec(name = 'ir', displayname = 'R')"
```

If you run the notebooks under a different kernel name, either select it in
Jupyter or override it on the command line:

```bash
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=ir \
    notebooks/01-getting-started.ipynb
```

The tutorials read data from a local `data/` directory and do not download
anything while running; see [notebooks/README.md](../notebooks/README.md) for the
exact data requirements and an example-data strategy.
`fdrtool` used to be required as an extra CRAN install; it is now pulled in transitively by other dependencies and you no longer need to install it explicitly.

## Verify the Installation

```bash
micromamba run -n contextnet R -e "
  library(ContextNet)
  check_ecosystem(verbose = TRUE)
"
```

This reports which optional packages (nichenetr, multinichenetr, Seurat, ccflowR) are available.

## Day-to-Day Usage

### Run R interactively

```bash
micromamba activate contextnet
R
```

### Run tests

```bash
# All tests
micromamba run -n contextnet R -e "devtools::test('.')"

# Specific test file (e.g., NicheNet tests)
micromamba run -n contextnet R -e "devtools::test('.', filter = 'nichenet')"

# Summary only (suppresses verbose output)
micromamba run -n contextnet R -e "devtools::test('.')" 2>&1 | tail -5
```

### Reinstall after code changes

```bash
# Fast reinstall (no dependency check)
micromamba run -n contextnet R CMD INSTALL --no-multiarch --with-keep.source .

# Full reinstall with dependency resolution
micromamba run -n contextnet R -e "remotes::install_local('.', dependencies = TRUE, upgrade = 'never')"
```

### Run package check

```bash
micromamba run -n contextnet R -e "devtools::check('.')"
```

### Launch Jupyter for tutorials

```bash
micromamba activate contextnet
jupyter notebook
# Open notebooks from the notebooks/ directory
```

## Managing the Environment

### List installed packages

```bash
micromamba list -n contextnet
```

### Search for a package

```bash
micromamba search -c conda-forge r-ggplot2
micromamba search -c bioconda bioconductor-limma
```

Package naming conventions:
- **CRAN packages:** `r-packagename` (lowercase), e.g., `r-dplyr`, `r-ggplot2`
- **Bioconductor packages:** `bioconductor-packagename`, e.g., `bioconductor-limma`
- **System libraries:** plain name, e.g., `nlopt`, `cairo`

### Install additional conda packages

```bash
micromamba install -n contextnet -c conda-forge r-somepackage -y
micromamba install -n contextnet -c bioconda bioconductor-somepackage -y
```

### Update packages

```bash
# Update a specific package
micromamba update -n contextnet -c conda-forge r-dplyr -y

# Update all packages (use with caution)
micromamba update -n contextnet --all -y
```

### Export the current environment

To capture exact versions for reproducibility:

```bash
micromamba list -n contextnet --export > environment-lock.txt
```

### Remove the environment

```bash
micromamba env remove -n contextnet -y
```

### Recreate from scratch

```bash
micromamba env remove -n contextnet -y
micromamba create -n contextnet -f environment.yml -y
# Then re-run the "Install Additional Dependencies" steps above
```

## Common Issues

| Problem | Solution |
|---------|----------|
| `BiocManager` not found when running step 2 manually | Already in `environment.yml` as `r-biocmanager`. Recreate the env or run `micromamba install -n contextnet -c conda-forge r-biocmanager -y`. |
| `hdf5r` configure error: "hdf5 does not seem to be installed" | Already in `environment.yml` as `r-hdf5r` (conda-prebuilt with bundled HDF5). Recreate the env or run `micromamba install -n contextnet -c conda-forge r-hdf5r -y`. |
| `remotes::install_github('settylab/ccflowR')` returns 404 | The repo is private. Export `GITHUB_PAT` (classic PAT with `repo` scope, or `$(gh auth token)`) and re-run the bootstrap. Without a PAT, ccflowR is skipped and only the tutorial Sankey visualization is unavailable. |
| `failed to lock directory` after a killed install | Stale `00LOCK-<pkg>` directories. Re-run `Rscript inst/bootstrap.R` — the bootstrap clears them. To clear manually: `rm -rf $CONDA_PREFIX/lib/R/library/00LOCK-*`. |
| `nloptr` fails to compile | Already in `environment.yml` as `nlopt + r-nloptr`. If your env predates that, run `micromamba install -n contextnet -c conda-forge nlopt r-nloptr r-lme4 -y`. |
| `repo.anaconda.com` unreachable | Use only `conda-forge` and `bioconda` channels; pass `--override-channels`. The shipped `environment.yml` already does this. |
| "Operation not permitted" on macOS source build | macOS extended attributes on the working tree: `xattr -cr /path/to/ContextNet`. |
| IRkernel not available in Jupyter | `micromamba run -n contextnet R -e "IRkernel::installspec()"`. |
| Package version conflicts | Pin the problematic package in `environment.yml` or install it last via R. |

## Quick Reference

```bash
# Create environment
micromamba create -n contextnet -f environment.yml -y

# Bootstrap GitHub-only deps and ContextNet (idempotent)
export GITHUB_PAT=$(gh auth token)   # required for private ccflowR
micromamba run -n contextnet Rscript inst/bootstrap.R

# Activate
micromamba activate contextnet

# Run R
micromamba run -n contextnet R

# Run tests
micromamba run -n contextnet R -e "devtools::test('.')"

# Verify
micromamba run -n contextnet R -e "library(ContextNet); check_ecosystem(verbose = TRUE)"

# Remove environment
micromamba env remove -n contextnet -y
```

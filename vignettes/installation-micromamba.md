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

Some packages are not available via conda and must be installed from R or GitHub.

### Activate the environment

```bash
micromamba activate contextnet
```

Or, to run commands without activating:

```bash
micromamba run -n contextnet R -e "<command>"
```

### Bioconductor packages (for MultiNicheNet)

```bash
micromamba run -n contextnet R -e "
  BiocManager::install(
    c('muscat', 'scater', 'sva', 'scran'),
    ask = FALSE, update = FALSE
  )
"
```

### GitHub packages

```bash
micromamba run -n contextnet R -e "
  remotes::install_github('saeyslab/nichenetr', upgrade = 'never')
  remotes::install_github('saeyslab/multinichenetr', upgrade = 'never')
  remotes::install_github('settylab/ccflowR', upgrade = 'never')
"
```

### CRAN packages not in conda

```bash
micromamba run -n contextnet R -e "install.packages('fdrtool')"
```

### Install ContextNet itself

```bash
micromamba run -n contextnet R -e "
  remotes::install_local('.', dependencies = TRUE, upgrade = 'never')
"
```

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
| `nloptr` fails to compile | `micromamba install -n contextnet -c conda-forge nlopt r-nloptr r-lme4 -y` |
| `fdrtool` not found in conda | Install from R: `install.packages("fdrtool")` |
| `repo.anaconda.com` unreachable | Use only `conda-forge` and `bioconda` channels; pass `--override-channels` |
| "Operation not permitted" on macOS | Clear extended attributes: `xattr -cr /path/to/ContextNet` |
| IRkernel not available in Jupyter | `micromamba run -n contextnet R -e "IRkernel::installspec()"` |
| Package version conflicts | Pin the problematic package in `environment.yml` or install it last via R |

## Quick Reference

```bash
# Create environment
micromamba create -n contextnet -f environment.yml -y

# Activate
micromamba activate contextnet

# Run R
micromamba run -n contextnet R

# Run tests
micromamba run -n contextnet R -e "devtools::test('.')"

# Install ContextNet
micromamba run -n contextnet R -e "remotes::install_local('.', dependencies = TRUE, upgrade = 'never')"

# Verify
micromamba run -n contextnet R -e "library(ContextNet); check_ecosystem(verbose = TRUE)"

# Remove environment
micromamba env remove -n contextnet -y
```

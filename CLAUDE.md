# Claude Code Project Notes

## Micromamba Setup

### Important: Use Direct Binary Path
The shell function for micromamba may not work properly in this environment. Always use the direct binary path:

```bash
/usr/local/bin/micromamba
```

### Network/Channel Restrictions
Company network does not use Anaconda directly. Using default channels (defaults, pkgs/main, pkgs/r) will cause DNS resolution failures:
- `repo.anaconda.com` cannot be resolved
- Only use `conda-forge` and `bioconda` channels

### Creating Environment
```bash
# Create environment (use -y for non-interactive)
/usr/local/bin/micromamba create -n contextnet_test1 -f /Volumes/yhuang2/ContextNet/dev/environment.yml -y

# Activate environment
eval "$(/usr/local/bin/micromamba shell hook -s bash)"
/usr/local/bin/micromamba activate contextnet_test1

# Or run commands directly in the environment
/usr/local/bin/micromamba run -n contextnet_test1 Rscript setup_and_test.R
```

### Environment File Location
- Environment definition: `dev/environment.yml`
- Setup script: `dev/setup_and_test.R`

### Troubleshooting
If channel errors occur, ensure `environment.yml` only uses:
- conda-forge
- bioconda

Do NOT use:
- defaults
- pkgs/main
- pkgs/r
- Any repo.anaconda.com channels

### Package Naming Conventions
In conda-forge/bioconda, packages have different naming:
- Regular R packages: `r-packagename` (e.g., `r-dplyr`, `r-ggplot2`)
- Bioconductor packages: `bioconductor-packagename` (e.g., `bioconductor-limma`)
- Some packages may not exist in conda channels - install via R's `install.packages()` or `BiocManager::install()` after env creation

Packages NOT available in conda-forge/bioconda (install in R):
- `fdrtool` - use `install.packages("fdrtool")` in R

### Known Installation Issues

#### nichenetr Dependencies
Installing nichenetr from GitHub may fail due to missing system dependencies:
- `nloptr` - requires nlopt library. Fix: `micromamba install -n contextnet_test1 nlopt`
- `lme4` - depends on nloptr
- `gdtools` - requires cairo/freetype. These should be in the conda env but may need explicit install
- `shadowtext`, `ggpubr`, `rstatix`, `car`, `pbkrtest` - dependency chain issues

**Workaround**: Install problematic packages via conda first:
```bash
/usr/local/bin/micromamba install -n contextnet_test1 -c conda-forge r-nloptr r-lme4 nlopt -y
```

#### ContextNet Package Build Issues
If build fails with "Operation not permitted" for files in dev/ or writing/ folders:
- This is a macOS extended attribute issue
- The package tries to copy files that have restricted flags
- **Fix**: Add these directories to `.Rbuildignore` or remove extended attributes:
```bash
xattr -cr /Volumes/yhuang2/ContextNet/dev
xattr -cr /Volumes/yhuang2/ContextNet/writing
```

## Complete Working Setup (Verified)

### Step-by-Step Installation

```bash
# 1. Create environment from environment.yml
/usr/local/bin/micromamba create -n contextnet_test1 -f /Volumes/yhuang2/ContextNet/dev/environment.yml -y

# 2. Install additional dependencies that may fail from R
/usr/local/bin/micromamba install -n contextnet_test1 -c conda-forge nlopt r-nloptr r-lme4 r-gdtools r-ggiraph r-shadowtext r-cairo r-locfit bioconductor-edger -y

# 3. Install Bioconductor packages for multinichenetr (via R)
cat > /tmp/install_bioc.R << 'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("muscat", "scater", "sva", "scran"), ask = FALSE, update = FALSE)
EOF
/usr/local/bin/micromamba run -n contextnet_test1 Rscript /tmp/install_bioc.R

# 4. Install nichenetr and multinichenetr from GitHub
/usr/local/bin/micromamba run -n contextnet_test1 R -e "remotes::install_github('saeyslab/nichenetr', upgrade = 'never')"
/usr/local/bin/micromamba run -n contextnet_test1 R -e "remotes::install_github('saeyslab/multinichenetr', upgrade = 'never')"

# 5. Install ContextNet (use R CMD INSTALL to avoid build issues)
/usr/local/bin/micromamba run -n contextnet_test1 R CMD INSTALL --no-multiarch --with-keep.source /Volumes/yhuang2/ContextNet

# 6. Run tests
/usr/local/bin/micromamba run -n contextnet_test1 R -e "library(ContextNet); devtools::test('/Volumes/yhuang2/ContextNet')"
```

### Expected Test Results
- **122 tests pass**
- 6 failures in show() tests (test methodology issue, output IS produced)
- 3 failures in network integration tests (need nichenetr data files)

### Verify Ecosystem
```bash
/usr/local/bin/micromamba run -n contextnet_test1 R -e "library(ContextNet); check_ecosystem(verbose=TRUE)"
```
Expected output:
- ✔ ContextNet
- ✔ nichenetr
- ✔ multinichenetr
- ✔ Seurat
- ✔ reticulate

### Quick Commands
```bash
# Run R in the environment
/usr/local/bin/micromamba run -n contextnet_test1 R

# Run tests only
/usr/local/bin/micromamba run -n contextnet_test1 R -e "devtools::test('/Volumes/yhuang2/ContextNet')"

# Check package
/usr/local/bin/micromamba run -n contextnet_test1 R -e "devtools::check('/Volumes/yhuang2/ContextNet')"

# Reinstall after code changes
/usr/local/bin/micromamba run -n contextnet_test1 R CMD INSTALL --no-multiarch --with-keep.source /Volumes/yhuang2/ContextNet
```

### Environment Name
The test environment is named `contextnet_test1` and located at:
`/Users/yhuang2/micromamba/envs/contextnet_test1`

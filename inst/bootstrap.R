# ContextNet bootstrap: installs GitHub-only and Bioconductor deps that
# `micromamba env create -f environment.yml` cannot supply, then installs
# ContextNet itself.
#
# Run AFTER creating and activating the conda env:
#
#   micromamba activate contextnet
#   Rscript inst/bootstrap.R
#
# Idempotent: skips packages that are already installed at a compatible
# version. Re-run safely after upgrading any single dependency.
#
# ccflowR (settylab/ccflowR) is a private repo. Export GITHUB_PAT before
# running to install it; without a PAT it is skipped with a NOTE and the
# Sankey-visualization paths in the tutorials are unavailable. Easiest
# way to mint one if you have `gh` logged in:
#
#   export GITHUB_PAT=$(gh auth token)
#   Rscript inst/bootstrap.R

# ── Config ────────────────────────────────────────────────────────────────────

# Pinned versions for upstream packages whose APIs drift between minor
# releases. `required = TRUE` means a missing/old version aborts the
# bootstrap; `required = FALSE` (e.g. ccflowR) just warns and continues
# so the install can succeed for users without access to the optional
# repo. Mirror these in DESCRIPTION's Suggests field.
#
# ccflowR is pinned by commit SHA, not version, because the package does
# not yet tag releases. The pin is captured at audit time and re-verified
# whenever the bootstrap is updated. See the install-audit report under
# `reports/` for the latest pin + reason.
PINS <- list(
    nichenetr      = list(repo = "saeyslab/nichenetr",      ref = NULL,
                          min  = "2.2.1", required = TRUE),
    multinichenetr = list(repo = "saeyslab/multinichenetr", ref = NULL,
                          min  = "2.1.0", required = TRUE),
    # `settylab/ccflowR` is a private repository. Users without read
    # access will get a 404 from `remotes::install_github`. The
    # bootstrap honours the GITHUB_PAT environment variable; without a
    # PAT, the install fails and ccflowR is skipped (every callsite is
    # namespace-guarded, so ContextNet itself works without it; only
    # the Sankey visualization in the tutorials is lost).
    ccflowR        = list(repo = "settylab/ccflowR",
                          ref  = "6056c0fe1d1cecfc68aa059bacfbe6bd16725231",
                          min  = NULL, required = FALSE)
)

# Bioconductor deps that multinichenetr requires but bioconda either does
# not ship or ships at incompatible versions.
BIOC_DEPS <- c("muscat", "scater", "sva", "scran")

# ── Helpers ───────────────────────────────────────────────────────────────────

# Local fallback for %||% so the script runs on R 4.1–4.3 (where base lacks
# it) without first having to install rlang.
`%||%` <- function(x, y) if (is.null(x)) y else x

msg <- function(...) cat(sprintf("[bootstrap] %s\n", paste0(...)))

has_min_version <- function(pkg, min_version) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(FALSE)
    if (is.null(min_version)) return(TRUE)
    utils::packageVersion(pkg) >= package_version(min_version)
}

ensure_installer <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        msg("installing ", pkg, " (build dependency)")
        utils::install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

# ── Stale lock cleanup ────────────────────────────────────────────────────────

# `install.packages()` leaves `00LOCK-<pkg>` directories behind if killed
# mid-build (e.g. SIGINT, OOM). On the next run the lock remains, every
# subsequent install of that package errors with "failed to lock
# directory", and the failure cascades to anything that Depends on it.
# Bootstrap is idempotent, so a re-run is the natural fix — but it has to
# clear the stale locks first, otherwise the re-run produces the same
# cascade. This is a no-op in the clean case.
.lib <- .libPaths()[1]
.locks <- list.files(.lib, pattern = "^00LOCK-", full.names = TRUE)
if (length(.locks)) {
    msg("clearing ", length(.locks), " stale 00LOCK director(y/ies) in ", .lib)
    unlink(.locks, recursive = TRUE, force = TRUE)
}

# ── Step 1: build dependencies ────────────────────────────────────────────────

msg("Step 1/4: ensure remotes + BiocManager are available")
ensure_installer("remotes")
ensure_installer("BiocManager")

# ── Step 2: Bioconductor deps for multinichenetr ──────────────────────────────

msg("Step 2/4: install Bioconductor deps for multinichenetr")
missing_bioc <- BIOC_DEPS[!vapply(BIOC_DEPS, requireNamespace,
                                  logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
    msg("  installing: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
    msg("  all Bioconductor deps already installed")
}

# ── Step 3: GitHub-only packages at pinned versions ───────────────────────────

msg("Step 3/4: install GitHub-only packages at pinned versions")
ccflowR_is_private <- !nzchar(Sys.getenv("GITHUB_PAT", unset = ""))
for (pkg in names(PINS)) {
    pin <- PINS[[pkg]]
    if (has_min_version(pkg, pin$min)) {
        installed <- as.character(utils::packageVersion(pkg))
        msg("  ", pkg, " ", installed, " already satisfies pin (>= ",
            pin$min %||% "any", ")")
        next
    }
    if (pkg == "ccflowR" && ccflowR_is_private) {
        msg("  NOTE: settylab/ccflowR is a private repository and ",
            "GITHUB_PAT is not set. Skipping ccflowR install. ",
            "To install: export GITHUB_PAT=<token-with-repo-scope> ",
            "and re-run this bootstrap. ContextNet works without ccflowR; ",
            "only the Sankey visualization in tutorial 02/03 will be unavailable.")
        next
    }
    spec <- if (!is.null(pin$ref)) paste0(pin$repo, "@", substr(pin$ref, 1, 7)) else pin$repo
    msg("  installing ", spec,
        if (!is.null(pin$min)) paste0(" (>= ", pin$min, ")") else "")
    install_ok <- tryCatch({
        remotes::install_github(pin$repo, ref = pin$ref %||% "HEAD",
                                upgrade = "never", quiet = FALSE)
        TRUE
    }, error = function(e) {
        msg("  WARN: install_github('", pin$repo, "') failed: ", conditionMessage(e))
        FALSE
    })
    if (!install_ok && pin$required) {
        stop(sprintf("[bootstrap] required package %s could not be installed from %s",
                     pkg, pin$repo))
    }
}

# Sanity-check that the post-install versions actually meet the pins.
# Optional pins (required = FALSE) are allowed to be missing.
for (pkg in names(PINS)) {
    pin <- PINS[[pkg]]
    if (has_min_version(pkg, pin$min)) next
    if (isTRUE(pin$required)) {
        stop(sprintf(
            "[bootstrap] %s installed but does not meet version pin (>= %s). Got %s.",
            pkg, pin$min %||% "any",
            tryCatch(as.character(utils::packageVersion(pkg)),
                     error = function(e) "<not installed>")
        ))
    } else {
        msg("  NOTE: optional package ", pkg, " not installed; skipping pin check")
    }
}

# ── Step 4: install ContextNet itself ─────────────────────────────────────────

# `inst/bootstrap.R` lives inside the package source tree. When run via
# `Rscript inst/bootstrap.R` the working directory is the package root.
msg("Step 4/4: install ContextNet from local source")
remotes::install_local(".", upgrade = "never", dependencies = TRUE,
                       quiet = FALSE)

# ── Verify ────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(ContextNet))
msg("ContextNet ", as.character(utils::packageVersion("ContextNet")), " installed")
ContextNet::check_ecosystem(verbose = TRUE)

msg("done. To run tests:")
msg("  Rscript -e 'devtools::test()'")

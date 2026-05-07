#' MultiNicheNet Interoperability Functions
#'
#' @description
#' Thin wrappers around multinichenetr functions for multi-condition
#' cell-cell communication analysis using context-specific networks.
#'
#' @name interop-multinichenet
NULL


#' Run MultiNicheNet multi-condition analysis
#'
#' @description A streamlined wrapper around multinichenetr's differential
#'   cell-cell communication workflow. Compares ligand activities between
#'   conditions using a ContextNetwork (with integrated GRN).
#'
#' @param object A Seurat object containing expression data with condition info
#' @param network A ContextNetwork object from integrate_networks()
#' @param sender_receiver_pairs List of sender-receiver pairs to analyze.
#'   Each element should be a character vector of length 2: c(sender, receiver)
#' @param condition_col Column name in metadata containing condition labels.
#'   REQUIRED.
#' @param contrast Named character vector specifying comparison:
#'   c(condition_oi = "Disease", condition_ref = "Control"). REQUIRED.
#' @param sample_col Column name in metadata containing sample/replicate IDs.
#'   REQUIRED for proper pseudobulk analysis.
#' @param group_col Column name in metadata containing cell type labels
#'   (default: "cell_type")
#' @param assay Assay to use from Seurat object (default: "RNA"). Used for
#'   expressed gene detection during ligand pre-filtering.
#' @param batches Column name(s) for batch correction, or NA for none (default: NA)
#' @param covariates Column name(s) for covariate correction, or NA for none (default: NA)
#' @param min_cells Minimum cells per sample/celltype/condition (default: 10)
#' @param logfc_threshold Log fold change threshold for DE (default: 0.5)
#' @param pval_threshold P-value threshold for DE (default: 0.05)
#' @param fraction_cutoff Minimum fraction of cells expressing a gene (default: 0.05)
#' @param top_n_target Top N target genes per ligand (default: 250)
#' @param empirical_pval Logical, use permutation-based empirical p-values for DE
#'   (default: FALSE). The multinichenetr wrapper defaults to TRUE, but this
#'   changes all p-values and produces different ligand rankings than the
#'   standard edgeR-based workflow shown in NicheNet/MultiNicheNet tutorials.
#'   Set to TRUE only if you specifically want empirical p-values.
#' @param lr_network Optional ligand-receptor network. If NULL, loads from
#'   nichenet_database or multinichenetr.
#' @param nichenet_database Path to directory containing NicheNet data files
#' @param n_cores Number of cores for parallel processing (default: 1)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item ligand_activities: data.frame of ligand activity scores per condition
#'     \item prioritization_table: data.frame of prioritized ligand-receptor pairs
#'     \item celltype_de: differential expression results per cell type
#'     \item multinichenet_output: raw output from multinichenetr (for advanced use)
#'   }
#'
#' @details
#' This function wraps multinichenetr's `multi_nichenet_analysis()` function,
#' which performs:
#' 1. Pseudobulk differential expression between conditions per cell type
#' 2. Ligand-receptor expression analysis
#' 3. Ligand activity prediction using your context-specific GRN
#' 4. Prioritization of differential ligand-receptor pairs
#'
#' The key advantage is that ligand-target predictions use your context-specific
#' GRN from integrate_networks().
#'
#' For more control over the analysis, use multinichenetr functions directly
#' with the ligand_target_matrix from your ContextNetwork object.
#'
#' @examples
#' \dontrun{
#' # Import and integrate GRN
#' grn <- import_grn("path/to/grn/")
#' network <- integrate_networks(grn, nichenet_database = "path/to/nichenet_data/")
#'
#' # Run MultiNicheNet analysis
#' result <- run_multinichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender_receiver_pairs = list(
#'         c("Macrophage", "Fibroblast"),
#'         c("T_cell", "Fibroblast")
#'     ),
#'     condition_col = "disease_status",
#'     contrast = c(condition_oi = "Disease", condition_ref = "Control"),
#'     sample_col = "sample_id"
#' )
#'
#' # Access results
#' head(result$prioritization_table)
#' }
#'
#' @export
run_multinichenet <- function(object,
                               network,
                               sender_receiver_pairs,
                               condition_col,
                               contrast,
                               sample_col,
                               group_col = "cell_type",
                               assay = "RNA",
                               batches = NA,
                               covariates = NA,
                               min_cells = 10,
                               logfc_threshold = 0.5,
                               pval_threshold = 0.05,
                               fraction_cutoff = 0.05,
                               top_n_target = 250,
                               empirical_pval = FALSE,
                               lr_network = NULL,
                               nichenet_database = NULL,
                               n_cores = 1,
                               verbose = TRUE) {
    # ── Check dependencies ──
    if (!requireNamespace("multinichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg multinichenetr} is required for run_multinichenet()",
            "i" = "Install with: {.code devtools::install_github('saeyslab/multinichenetr')}"
        ))
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg Seurat} is required for run_multinichenet()",
            "i" = "Install with: {.code install.packages('Seurat')}"
        ))
    }

    # ── Validate inputs ──
    # Validate network first (doesn't require Seurat)
    if (!inherits(network, "ContextNetwork")) {
        cli::cli_abort(c(
            "{.arg network} must be a ContextNetwork object",
            "i" = "Create one with {.fn integrate_networks}"
        ))
    }

    # Validate required parameters before Seurat check
    # This allows testing parameter validation without a Seurat object
    if (is.null(sender_receiver_pairs)) {
        cli::cli_abort(c(
            "{.arg sender_receiver_pairs} is required",
            "i" = "Provide a list of c(sender, receiver) pairs"
        ))
    }
    .validate_sender_receiver_pairs(sender_receiver_pairs)

    if (is.null(condition_col) || !nzchar(condition_col)) {
        cli::cli_abort(c(
            "{.arg condition_col} is required",
            "i" = "Specify the metadata column containing condition labels"
        ))
    }

    if (is.null(contrast)) {
        cli::cli_abort(c(
            "{.arg contrast} is required",
            "i" = "Provide: c(condition_oi = 'Disease', condition_ref = 'Control')"
        ))
    }
    .validate_contrast(contrast)

    if (missing(sample_col) || is.null(sample_col)) {
        cli::cli_abort(c(
            "{.arg sample_col} is required for multinichenetr analysis",
            "i" = "Specify the metadata column containing sample/replicate IDs"
        ))
    }

    # Now validate Seurat object
    if (!inherits(object, "Seurat")) {
        cli::cli_abort(c(
            "{.arg object} must be a Seurat object"
        ))
    }

    # Validate columns exist in metadata
    required_cols <- c(group_col, condition_col, sample_col)
    # Add batch columns if specified
    if (!all(is.na(batches))) {
        required_cols <- c(required_cols, batches[!is.na(batches)])
    }
    # Add covariate columns if specified
    if (!all(is.na(covariates))) {
        required_cols <- c(required_cols, covariates[!is.na(covariates)])
    }

    missing_cols <- setdiff(required_cols, colnames(object@meta.data))
    if (length(missing_cols) > 0) {
        cli::cli_abort(c(
            "Required columns not found in Seurat metadata",
            "x" = "Missing: {.val {missing_cols}}",
            "i" = "Available: {.val {colnames(object@meta.data)}}"
        ))
    }

    # Validate cell types exist
    available_celltypes <- unique(object@meta.data[[group_col]])
    all_celltypes <- unique(unlist(sender_receiver_pairs))
    missing_celltypes <- setdiff(all_celltypes, available_celltypes)
    if (length(missing_celltypes) > 0) {
        cli::cli_abort(c(
            "Cell types not found in metadata",
            "x" = "Missing: {.val {missing_celltypes}}",
            "i" = "Available: {.val {available_celltypes}}"
        ))
    }

    # Validate conditions exist
    available_conditions <- unique(object@meta.data[[condition_col]])
    required_conditions <- c(contrast["condition_oi"], contrast["condition_ref"])
    missing_conditions <- setdiff(required_conditions, available_conditions)
    if (length(missing_conditions) > 0) {
        cli::cli_abort(c(
            "Conditions not found in metadata",
            "x" = "Missing: {.val {missing_conditions}}",
            "i" = "Available: {.val {available_conditions}}"
        ))
    }

    # ── Load lr_network if needed ──
    if (is.null(lr_network)) {
        # Try 1: Load from local nichenet_database (prefer allInfo version)
        if (!is.null(nichenet_database)) {
            lr_files <- c(
                "lr_network_human_allInfo_30112033.rds",
                "lr_network_human_21122021.rds",
                "lr_network_human.rds"
            )
            for (f in lr_files) {
                path <- file.path(nichenet_database, f)
                if (file.exists(path)) {
                    lr_network <- readRDS(path)
                    if (verbose) cli::cli_alert_success("Loaded lr_network from {.file {f}}")
                    break
                }
            }
        }

        # Try 2: Download newest lr_network from Zenodo
        if (is.null(lr_network)) {
            lr_network <- tryCatch({
                if (verbose) cli::cli_alert("Downloading lr_network from Zenodo...")
                lr_all <- readRDS(url(
                    "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
                ))
                if (verbose) cli::cli_alert_success("Downloaded lr_network from Zenodo")
                # Save to nichenet_database for future use if path is available
                if (!is.null(nichenet_database) && dir.exists(nichenet_database)) {
                    save_path <- file.path(nichenet_database,
                                            "lr_network_human_allInfo_30112033.rds")
                    tryCatch({
                        saveRDS(lr_all, save_path)
                        if (verbose) cli::cli_alert_info("Cached lr_network to {.file {save_path}}")
                    }, error = function(e) NULL)
                }
                lr_all
            }, error = function(e) {
                if (verbose) cli::cli_alert_warning(
                    "Zenodo download failed: {e$message}"
                )
                NULL
            })
        }

        # Try 3: Fall back to package data
        if (is.null(lr_network)) {
            lr_network <- tryCatch({
                multinichenetr::lr_network_human
            }, error = function(e) {
                tryCatch({
                    nichenetr::lr_network
                }, error = function(e2) {
                    cli::cli_abort(c(
                        "Cannot load lr_network",
                        "i" = "Provide via {.arg lr_network} or {.arg nichenet_database}",
                        "i" = "Or use {.fn load_nichenet_data} to download from Zenodo"
                    ))
                })
            })
        }
    }

    # Ensure lr_network has correct column names (ligand, receptor)
    if ("from" %in% colnames(lr_network) && !"ligand" %in% colnames(lr_network)) {
        lr_network <- lr_network %>%
            dplyr::rename(ligand = "from", receptor = "to")
    }

    # Convert gene aliases to official HGNC symbols in lr_network.
    # The Zenodo lr_network and SCE gene names both benefit from alias resolution
    # (e.g., YARS → YARS1). Without this, some LR pairs won't match SCE genes.
    tryCatch({
        lr_network <- lr_network %>%
            dplyr::mutate(
                ligand = nichenetr::convert_alias_to_symbols(ligand, organism = "human"),
                receptor = nichenetr::convert_alias_to_symbols(receptor, organism = "human")
            )
    }, error = function(e) {
        if (verbose) {
            cli::cli_alert_warning(
                "Could not run convert_alias_to_symbols on lr_network: {e$message}"
            )
        }
    })

    # Normalize lr_network gene names to match SCE conventions.
    # alias_to_symbol_SCE() + makenames_SCE() converts SCE gene names to
    # official HGNC symbols with valid R names (e.g., "HLA-A" → "HLA.A").
    # The lr_network must undergo the same transformation so LR pairs match.
    lr_network <- lr_network %>%
        dplyr::mutate(
            ligand = make.names(ligand),
            receptor = make.names(receptor)
        ) %>%
        dplyr::distinct(ligand, receptor)

    # ── Extract ligand-target matrix from network ──
    ligand_target_matrix <- network@ligand_target_matrix

    # ── Filter LTM against lr_network (ensures consistency) ──
    # Keep only ligands that appear in both LTM and lr_network. This matches
    # the manual multinichenetr workflow which does:
    #   lr_network = lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
    #   ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]
    # Ligands not in lr_network have no known receptor, so they cannot participate
    # in ligand-receptor interactions. This is NOT the same as pre-filtering by
    # expression (which was reverted) — it's just removing impossible interactions.
    lr_ligands <- unique(lr_network$ligand)
    ltm_ligands <- colnames(ligand_target_matrix)
    shared_ligands <- intersect(lr_ligands, ltm_ligands)

    lr_network <- lr_network %>%
        dplyr::filter(ligand %in% shared_ligands)
    ligand_target_matrix <- ligand_target_matrix[, shared_ligands, drop = FALSE]

    if (verbose) {
        cli::cli_alert_info(
            "LTM filtered: {.val {length(shared_ligands)}} ligands shared with lr_network (from {.val {length(ltm_ligands)}})"
        )
    }

    if (verbose) {
        cli::cli_h2("MultiNicheNet Differential Analysis")
        cli::cli_alert_info("Condition of interest: {.val {contrast['condition_oi']}}")
        cli::cli_alert_info("Reference condition: {.val {contrast['condition_ref']}}")
        cli::cli_alert_info("Sender-receiver pairs: {.val {length(sender_receiver_pairs)}}")
    }

    # ── Prepare sender-receiver info ──
    senders <- unique(sapply(sender_receiver_pairs, `[`, 1))
    receivers <- unique(sapply(sender_receiver_pairs, `[`, 2))

    if (verbose) {
        cli::cli_alert_info("Senders: {.val {senders}}")
        cli::cli_alert_info("Receivers: {.val {receivers}}")
    }

    # ── Subset to relevant cell types (performance optimization) ──
    # Mirrors run_nichenet approach: only process cells from sender/receiver types.
    # This avoids computing pseudobulk DE for irrelevant cell types.
    relevant_celltypes <- unique(c(senders, receivers))
    cells_keep <- object@meta.data[[group_col]] %in% relevant_celltypes
    n_total_cells <- ncol(object)
    object <- object[, cells_keep]

    if (verbose) {
        cli::cli_alert_info(
            "Subset to sender/receiver cells: {.val {ncol(object)}} / {.val {n_total_cells}}"
        )
    }

    # ── Pre-filter cell types with insufficient cells for pseudobulk ──
    # multinichenetr requires >= 2 samples with >= min_cells per condition.
    # Cell types that fail this are excluded by pbDS with an error. Pre-filter
    # them here to avoid wasting time and producing confusing error messages.
    condition_oi_val <- contrast["condition_oi"]
    condition_ref_val <- contrast["condition_ref"]
    meta <- object@meta.data
    insufficient_celltypes <- character()

    for (ct in relevant_celltypes) {
        ct_meta <- meta[meta[[group_col]] == ct, ]
        for (cond in c(condition_oi_val, condition_ref_val)) {
            ct_cond <- ct_meta[ct_meta[[condition_col]] == cond, ]
            samples_with_enough <- sum(
                table(ct_cond[[sample_col]]) >= min_cells
            )
            if (samples_with_enough < 2) {
                insufficient_celltypes <- c(insufficient_celltypes, ct)
                break
            }
        }
    }

    if (length(insufficient_celltypes) > 0) {
        insufficient_celltypes <- unique(insufficient_celltypes)
        cli::cli_warn(c(
            "Excluding cell types with insufficient cells for pseudobulk DE",
            "x" = "Excluded: {.val {insufficient_celltypes}}",
            "i" = "Requires >= 2 samples with >= {min_cells} cells per condition"
        ))
        cells_keep2 <- !(object@meta.data[[group_col]] %in% insufficient_celltypes)
        object <- object[, cells_keep2]
        senders <- setdiff(senders, insufficient_celltypes)
        receivers <- setdiff(receivers, insufficient_celltypes)

        if (length(senders) == 0) {
            cli::cli_abort("No sender cell types remain after filtering")
        }
        if (length(receivers) == 0) {
            cli::cli_abort("No receiver cell types remain after filtering")
        }
    }

    # NOTE: Unlike run_nichenet(), we do NOT pre-filter the ligand-target matrix
    # by expression here (e.g., removing ligands not expressed in sender cells).
    # multinichenetr computes ligand activities for ALL cell types as receivers
    # (not just specified ones), and its prioritization_score includes sub-scores
    # normalized relative to all tested ligands (e.g., activity_scaled). Filtering
    # by expression would change these normalizations and produce different rankings.
    # The LTM IS filtered against lr_network above (removing ligands with no known
    # receptor), which is a structural filter matching the manual workflow.

    # ── Convert contrast to multinichenetr format (bidirectional) ──
    # multinichenetr expects both directions for proper activity_scaled normalization
    # and prioritization_score computation.
    # Our format: c(condition_oi = "Disease", condition_ref = "Control")
    # multinichenetr format: contrasts_oi = "'Disease-Control','Control-Disease'"
    #                        contrast_tbl with both directions
    condition_oi <- contrast["condition_oi"]
    condition_ref <- contrast["condition_ref"]
    contrasts_oi <- paste0(
        "'", condition_oi, "-", condition_ref, "',",
        "'", condition_ref, "-", condition_oi, "'"
    )
    contrast_tbl <- data.frame(
        contrast = c(
            paste0(condition_oi, "-", condition_ref),
            paste0(condition_ref, "-", condition_oi)
        ),
        group = c(condition_oi, condition_ref),
        stringsAsFactors = FALSE
    )

    # ── Sanitize cell type names (multinichenetr requires valid R names) ──
    original_celltypes <- unique(as.character(object@meta.data[[group_col]]))
    safe_celltypes <- make.names(original_celltypes)
    needs_renaming <- any(original_celltypes != safe_celltypes)
    reverse_map <- NULL

    if (needs_renaming) {
        name_map <- setNames(safe_celltypes, original_celltypes)
        reverse_map <- setNames(original_celltypes, safe_celltypes)
        if (verbose) {
            changed <- original_celltypes[original_celltypes != safe_celltypes]
            cli::cli_alert_info(
                "Renaming cell types for multinichenetr compatibility: {.val {changed}}"
            )
        }
        object@meta.data[[group_col]] <- name_map[as.character(object@meta.data[[group_col]])]
        senders <- name_map[senders]
        receivers <- name_map[receivers]
    }

    # ── Convert Seurat to SingleCellExperiment ──
    sce <- Seurat::as.SingleCellExperiment(object)

    # Normalize gene names to match NicheNet/LTM conventions:
    # - alias_to_symbol_SCE: maps gene aliases to official HGNC symbols
    # - makenames_SCE: converts to valid R names (e.g., "NKX3-1" -> "NKX3.1")
    # Without this, gene name mismatches between SCE and LTM reduce the
    # available gene set for DE and ligand activity computation.
    tryCatch({
        sce <- multinichenetr::alias_to_symbol_SCE(sce, "human")
    }, error = function(e) {
        if (verbose) {
            cli::cli_warn(c(
                "Could not run alias_to_symbol_SCE: {e$message}",
                "i" = "Gene alias conversion skipped; some gene names may not match LTM"
            ))
        }
    })
    sce <- multinichenetr::makenames_SCE(sce)

    # ── Run multinichenetr analysis ──
    if (verbose) cli::cli_alert("Running multi_nichenet_analysis()...")

    multinichenet_output <- multinichenetr::multi_nichenet_analysis(
        sce = sce,
        celltype_id = group_col,
        sample_id = sample_col,
        group_id = condition_col,
        batches = batches,
        covariates = covariates,
        lr_network = lr_network,
        ligand_target_matrix = ligand_target_matrix,
        contrasts_oi = contrasts_oi,
        contrast_tbl = contrast_tbl,
        senders_oi = senders,
        receivers_oi = receivers,
        fraction_cutoff = fraction_cutoff,
        min_cells = min_cells,
        logFC_threshold = logfc_threshold,
        p_val_threshold = pval_threshold,
        empirical_pval = empirical_pval,
        top_n_target = top_n_target,
        verbose = verbose,
        n.cores = n_cores
    )

    # ── Validate upstream output structure ──
    # Catches silent API changes in multinichenetr (e.g., renamed/removed fields)
    .assert_upstream_output(
        multinichenet_output,
        required_elements = c(
            "ligand_activities_targets_DEgenes",
            "prioritization_tables",
            "celltype_de"
        ),
        caller = "multinichenetr::multi_nichenet_analysis"
    )

    # Validate nested structure — core fields are hard errors
    .assert_upstream_output(
        multinichenet_output$ligand_activities_targets_DEgenes,
        required_elements = "ligand_activities",
        caller = "multi_nichenet_analysis()$ligand_activities_targets_DEgenes"
    )

    # Prioritization table: multinichenetr >= 2.1.0 renamed
    # "prioritization_tbl_ligand_receptor" → "group_prioritization_tbl".
    # Try new name first, fall back to old name.
    prio_tables <- multinichenet_output$prioritization_tables
    prio_tbl <- prio_tables$group_prioritization_tbl %||%
        prio_tables$prioritization_tbl_ligand_receptor

    if (!is.null(prio_tables) && is.null(prio_tbl)) {
        cli::cli_warn(c(
            "multinichenetr output missing prioritization table",
            "i" = "Expected {.field group_prioritization_tbl} or {.field prioritization_tbl_ligand_receptor}",
            "i" = "Got elements: {.val {names(prio_tables)}}",
            "i" = "Ligand ranking will fall back to ligand_activities"
        ))
    }

    # ── Prepare output ──
    result <- list(
        ligand_activities = multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities,
        prioritization_table = prio_tbl,
        celltype_de = multinichenet_output$celltype_de,
        multinichenet_output = multinichenet_output
    )

    # Restore original cell type names if they were sanitized
    if (needs_renaming) {
        result <- .restore_celltype_names(result, reverse_map)
    }

    if (verbose) {
        cli::cli_alert_success("Analysis complete")
        if (!is.null(result$prioritization_table) && nrow(result$prioritization_table) > 0) {
            cli::cli_dl(c(
                "Prioritized LR pairs" = "{.val {nrow(result$prioritization_table)}}",
                "Unique ligands" = "{.val {length(unique(result$prioritization_table$ligand))}}"
            ))
        }
    }

    result
}


#' Validate sender-receiver pairs structure
#'
#' @param pairs List of sender-receiver pairs
#'
#' @keywords internal
.validate_sender_receiver_pairs <- function(pairs) {
    if (!is.list(pairs)) {
        cli::cli_abort("{.arg sender_receiver_pairs} must be a list")
    }

    for (i in seq_along(pairs)) {
        pair <- pairs[[i]]
        if (!is.character(pair) || length(pair) != 2) {
            cli::cli_abort(c(
                "Each sender-receiver pair must be a character vector of length 2",
                "x" = "Element {i} has length {length(pair)}"
            ))
        }
    }

    invisible(TRUE)
}


#' Validate contrast structure
#'
#' @param contrast Named character vector
#'
#' @keywords internal
.validate_contrast <- function(contrast) {
    required <- c("condition_oi", "condition_ref")
    missing <- setdiff(required, names(contrast))

    if (length(missing) > 0) {
        cli::cli_abort(c(
            "{.arg contrast} missing required elements",
            "x" = "Missing: {.val {missing}}",
            "i" = "Expected: c(condition_oi = 'Disease', condition_ref = 'Control')"
        ))
    }

    invisible(TRUE)
}


#' Restore original cell type names in multinichenet output
#'
#' @description multinichenetr requires syntactically valid R names for cell type
#'   factors. This function recursively maps sanitized names (e.g., "T.Cell")
#'   back to their originals (e.g., "T Cell") in the result structure.
#'
#' @param x A list, data.frame, or other object from multinichenet output
#' @param reverse_map Named character vector: sanitized names -> original names
#'
#' @return The input with cell type names restored
#' @keywords internal
.restore_celltype_names <- function(x, reverse_map) {
    if (is.data.frame(x)) {
        for (col in colnames(x)) {
            if (is.character(x[[col]])) {
                mapped <- reverse_map[x[[col]]]
                has_match <- !is.na(mapped)
                if (any(has_match)) {
                    x[[col]][has_match] <- mapped[has_match]
                }
            } else if (is.factor(x[[col]])) {
                lvls <- levels(x[[col]])
                mapped <- reverse_map[lvls]
                has_match <- !is.na(mapped)
                if (any(has_match)) {
                    lvls[has_match] <- mapped[has_match]
                    levels(x[[col]]) <- lvls
                }
            }
        }
    } else if (is.list(x) && !isS4(x) && !is.environment(x)) {
        # Map list names
        if (!is.null(names(x))) {
            nm <- names(x)
            mapped <- reverse_map[nm]
            has_match <- !is.na(mapped)
            if (any(has_match)) {
                nm[has_match] <- mapped[has_match]
                names(x) <- nm
            }
        }
        # Recurse into list elements, skipping elements that can't be accessed
        # (multinichenet output may contain model objects, formulas, etc.)
        for (i in seq_along(x)) {
            tryCatch({
                x[[i]] <- .restore_celltype_names(x[[i]], reverse_map)
            }, error = function(e) NULL)
        }
    }
    x
}

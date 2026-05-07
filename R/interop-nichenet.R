#' NicheNet Interoperability Functions
#'
#' @description
#' Thin wrappers around nichenetr functions for streamlined cell-cell
#' communication analysis using context-specific networks.
#'
#' @name interop-nichenet
NULL


#' Run NicheNet ligand activity analysis
#'
#' @description A streamlined wrapper around nichenetr's ligand activity
#'   prediction workflow. Uses a ContextNetwork (with integrated GRN) to
#'   predict which ligands from sender cells regulate target genes in
#'   receiver cells.
#'
#' @param object A Seurat object containing expression data
#' @param network A ContextNetwork object from integrate_networks()
#' @param sender Character, name of sender cell type(s) in metadata
#' @param receiver Character, name of receiver cell type in metadata
#' @param geneset_oi Character vector of genes of interest in receiver
#'   (e.g., differentially expressed genes). REQUIRED.
#' @param group_col Column name in metadata containing cell type labels
#'   (default: "cell_type")
#' @param assay Assay to use from Seurat object (default: "RNA")
#' @param slot Slot to use for expression data (default: "data")
#' @param min_pct Minimum percentage of cells expressing a gene to consider
#'   it expressed (default: 0.10)
#' @param top_n_ligands Number of top ligands to return (default: 20)
#' @param ligand_target_cutoff Minimum regulatory potential score to include a
#'   ligand-target link (default: 0)
#' @param lr_network Optional ligand-receptor network. If NULL, uses
#'   nichenetr's default or loads from nichenet_database.
#' @param nichenet_database Path to directory containing NicheNet data files
#'   (lr_network_human.rds, etc.). If NULL, uses nichenetr package data.
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item ligand_activities: data.frame of ligand activity scores
#'     \item ligand_target_links: data.frame of ligand-target regulatory potential
#'     \item top_ligands: character vector of top-ranked ligands
#'     \item receptor_info: data.frame of ligand-receptor pairs
#'     \item nichenet_input: list of inputs used (for reproducibility)
#'   }
#'
#' @details
#' This function performs the standard NicheNet ligand activity analysis:
#' 1. Identifies genes expressed in sender and receiver cells
#' 2. Finds potential ligands (expressed in sender, receptor in receiver)
#' 3. Predicts ligand activities based on geneset_oi
#' 4. Extracts ligand-target regulatory potential
#'
#' The key difference from standard NicheNet is that the ligand-target
#' matrix comes from integrate_networks(), which incorporates your
#' context-specific GRN.
#'
#' @examples
#' \dontrun{
#' # Import and integrate GRN
#' grn <- import_grn("path/to/grn/")
#' network <- integrate_networks(grn, nichenet_database = "path/to/nichenet_data/")
#'
#' # Run NicheNet analysis
#' result <- run_nichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = "Macrophage",
#'     receiver = "Fibroblast",
#'     geneset_oi = de_genes,
#'     group_col = "cell_type"
#' )
#'
#' # Access results
#' head(result$ligand_activities)
#' result$top_ligands
#' }
#'
#' @export
run_nichenet <- function(object,
                          network,
                          sender,
                          receiver,
                          geneset_oi,
                          group_col = "cell_type",
                          assay = "RNA",
                          slot = "data",
                          min_pct = 0.10,
                          top_n_ligands = 20,
                          ligand_target_cutoff = 0,
                          lr_network = NULL,
                          nichenet_database = NULL,
                          verbose = TRUE) {
    # ── Check dependencies ──
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg nichenetr} is required for run_nichenet()",
            "i" = "Install with: {.code devtools::install_github('saeyslab/nichenetr')}"
        ))
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg Seurat} is required for run_nichenet()",
            "i" = "Install with: {.code install.packages('Seurat')}"
        ))
    }

    # ── Validate inputs ──
    if (!inherits(network, "ContextNetwork")) {
        cli::cli_abort(c(
            "{.arg network} must be a ContextNetwork object",
            "i" = "Create one with {.fn integrate_networks}"
        ))
    }

    if (!inherits(object, "Seurat")) {
        cli::cli_abort(c(
            "{.arg object} must be a Seurat object"
        ))
    }

    if (is.null(geneset_oi) || length(geneset_oi) == 0) {
        cli::cli_abort(c(
            "{.arg geneset_oi} is required",
            "i" = "Provide a character vector of genes of interest (e.g., DE genes)"
        ))
    }

    # Validate cell types exist in metadata
    if (!group_col %in% colnames(object@meta.data)) {
        cli::cli_abort(c(
            "Column {.val {group_col}} not found in Seurat metadata",
            "i" = "Available columns: {.val {colnames(object@meta.data)}}"
        ))
    }

    available_celltypes <- unique(object@meta.data[[group_col]])
    missing_sender <- setdiff(sender, available_celltypes)
    if (length(missing_sender) > 0) {
        cli::cli_abort(c(
            "Sender cell type(s) not found in metadata",
            "x" = "Missing: {.val {missing_sender}}",
            "i" = "Available: {.val {available_celltypes}}"
        ))
    }

    if (!receiver %in% available_celltypes) {
        cli::cli_abort(c(
            "Receiver cell type not found in metadata",
            "x" = "Missing: {.val {receiver}}",
            "i" = "Available: {.val {available_celltypes}}"
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
                nichenetr::lr_network
            }, error = function(e) {
                cli::cli_abort(c(
                    "Cannot load lr_network",
                    "i" = "Provide via {.arg lr_network} or {.arg nichenet_database}",
                    "i" = "Or use {.fn load_nichenet_data} to download from Zenodo"
                ))
            })
        }
    }

    # Ensure lr_network has correct column names (ligand, receptor)
    if ("from" %in% colnames(lr_network) && !"ligand" %in% colnames(lr_network)) {
        lr_network <- lr_network %>%
            dplyr::rename(ligand = "from", receptor = "to")
    }

    # Convert gene aliases to official HGNC symbols in lr_network
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

    # Normalize gene names to valid R names and deduplicate
    lr_network <- lr_network %>%
        dplyr::mutate(
            ligand = make.names(ligand),
            receptor = make.names(receptor)
        ) %>%
        dplyr::distinct(ligand, receptor)

    # ── Extract ligand-target matrix from network ──
    ligand_target_matrix <- network@ligand_target_matrix

    # ── Filter LTM against lr_network ──
    # Keep only ligands that appear in both LTM and lr_network. Ligands not in
    # lr_network have no known receptor and cannot participate in L-R interactions.
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
        cli::cli_h2("NicheNet Ligand Activity Analysis")
        cli::cli_alert_info("Sender: {.val {paste(sender, collapse = ', ')}}")
        cli::cli_alert_info("Receiver: {.val {receiver}}")
        cli::cli_alert_info("Genes of interest: {.val {length(geneset_oi)}}")
    }

    # ── Extract expression data ──
    if (verbose) cli::cli_alert("Extracting expression data...")

    expression_data <- .extract_expression_nichenet(
        object = object,
        sender = sender,
        receiver = receiver,
        group_col = group_col,
        assay = assay,
        slot = slot,
        min_pct = min_pct
    )

    sender_expressed <- expression_data$sender_expressed
    receiver_expressed <- expression_data$receiver_expressed
    background_genes <- expression_data$background_genes

    if (verbose) {
        cli::cli_alert_info("Sender expressed genes: {.val {length(sender_expressed)}}")
        cli::cli_alert_info("Receiver expressed genes: {.val {length(receiver_expressed)}}")
    }

    # ── Get potential ligands ──
    potential_ligands <- .get_potential_ligands(
        lr_network = lr_network,
        sender_expressed = sender_expressed,
        receiver_expressed = receiver_expressed
    )

    if (length(potential_ligands) == 0) {
        cli::cli_abort(c(
            "No potential ligands found",
            "i" = "Check that sender expresses ligands with receptors in receiver"
        ))
    }

    if (verbose) {
        cli::cli_alert_info("Potential ligands: {.val {length(potential_ligands)}}")
    }

    # ── Filter ligands to those in matrix ──
    ligands_in_matrix <- intersect(potential_ligands, colnames(ligand_target_matrix))

    if (length(ligands_in_matrix) == 0) {
        cli::cli_abort(c(
            "None of the potential ligands are in the ligand-target matrix",
            "i" = "This may indicate a mismatch between lr_network and the integrated network"
        ))
    }

    if (verbose && length(ligands_in_matrix) < length(potential_ligands)) {
        cli::cli_alert_warning(
            "Only {length(ligands_in_matrix)}/{length(potential_ligands)} ligands in matrix"
        )
    }

    # ── Filter geneset_oi to genes in matrix ──
    geneset_oi_filtered <- intersect(geneset_oi, rownames(ligand_target_matrix))

    if (length(geneset_oi_filtered) == 0) {
        cli::cli_abort(c(
            "None of the genes of interest are in the ligand-target matrix",
            "i" = "Check gene naming (e.g., symbols vs ENSEMBL IDs)"
        ))
    }

    if (verbose && length(geneset_oi_filtered) < length(geneset_oi)) {
        cli::cli_alert_warning(
            "Only {length(geneset_oi_filtered)}/{length(geneset_oi)} geneset genes in matrix"
        )
    }

    # ── Predict ligand activities ──
    if (verbose) cli::cli_alert("Predicting ligand activities...")

    # Background genes: expressed in receiver, in matrix
    background_filtered <- intersect(background_genes, rownames(ligand_target_matrix))

    ligand_activities <- nichenetr::predict_ligand_activities(
        geneset = geneset_oi_filtered,
        background_expressed_genes = background_filtered,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = ligands_in_matrix
    )

    # Validate upstream output — catches silent API changes
    .assert_upstream_output(
        ligand_activities,
        required_cols = c("test_ligand", "aupr_corrected"),
        caller = "nichenetr::predict_ligand_activities"
    )

    # Sort by activity score
    ligand_activities <- ligand_activities[order(
        ligand_activities$aupr_corrected,
        decreasing = TRUE
    ), ]

    # ── Get top ligands ──
    top_ligands <- utils::head(ligand_activities$test_ligand, top_n_ligands)

    if (verbose) {
        cli::cli_alert_success("Top {length(top_ligands)} ligands identified")
    }

    # ── Extract ligand-target links ──
    if (verbose) cli::cli_alert("Extracting ligand-target regulatory potential...")

    # Get active ligand-target links using get_active_ligand_target_df
    # (replaces deprecated get_active_ligand_target_links)
    # Subset ligand_target_matrix to top ligands first
    top_ligands_in_matrix <- intersect(top_ligands, colnames(ligand_target_matrix))
    if (length(top_ligands_in_matrix) > 0) {
        ltm_subset <- ligand_target_matrix[, top_ligands_in_matrix, drop = FALSE]

        # Create named logical vector for response
        # Names are all genes (rows of ltm), TRUE for genes in geneset_oi
        all_targets <- rownames(ltm_subset)
        response_vec <- stats::setNames(
            all_targets %in% geneset_oi_filtered,
            all_targets
        )

        active_ligand_target_links <- nichenetr::get_active_ligand_target_df(
            response = response_vec,
            ligand_target_matrix = ltm_subset,
            ligands_position = "cols",
            cutoff = ligand_target_cutoff
        )

        # Validate upstream output — catches column renames across versions
        .assert_upstream_output(
            active_ligand_target_links,
            required_cols = c("ligand", "target"),
            caller = "nichenetr::get_active_ligand_target_df"
        )

        # Rename 'score' to 'weight' for compatibility with nichenetr plotting
        # (prepare_ligand_target_visualization expects 'weight')
        colnames(active_ligand_target_links)[colnames(active_ligand_target_links) == "score"] <- "weight"
    } else {
        active_ligand_target_links <- data.frame(
            ligand = character(),
            target = character(),
            weight = numeric(),
            stringsAsFactors = FALSE
        )
    }

    # ── Get receptor information ──
    receptor_info <- lr_network[lr_network$ligand %in% top_ligands, ]
    receptor_info <- receptor_info[receptor_info$receptor %in% receiver_expressed, ]
    # Rename to from/to for ContextResult compatibility
    colnames(receptor_info)[colnames(receptor_info) == "ligand"] <- "from"
    colnames(receptor_info)[colnames(receptor_info) == "receptor"] <- "to"

    # ── Prepare output ──
    result <- list(
        ligand_activities = ligand_activities,
        ligand_target_links = active_ligand_target_links,
        top_ligands = top_ligands,
        receptor_info = receptor_info,
        nichenet_input = list(
            geneset_oi = geneset_oi_filtered,
            background_genes = background_filtered,
            potential_ligands = ligands_in_matrix,
            sender = sender,
            receiver = receiver,
            sender_expressed = sender_expressed,
            receiver_expressed = receiver_expressed
        )
    )

    if (verbose) {
        cli::cli_alert_success("Analysis complete")
        cli::cli_dl(c(
            "Top ligand" = "{.val {top_ligands[1]}} (AUPR: {round(ligand_activities$aupr_corrected[1], 4)})",
            "Ligand-target links" = "{.val {nrow(active_ligand_target_links)}}"
        ))
    }

    result
}


#' Extract expression data from Seurat object
#'
#' @param object Seurat object
#' @param sender Sender cell type(s)
#' @param receiver Receiver cell type
#' @param group_col Column with cell type labels
#' @param assay Assay to use
#' @param slot Slot to use
#' @param min_pct Minimum percent expression
#'
#' @return List with sender_expressed, receiver_expressed, background_genes
#'
#' @keywords internal
.extract_expression_nichenet <- function(object,
                                          sender,
                                          receiver,
                                          group_col,
                                          assay,
                                          slot,
                                          min_pct) {
    # Use nichenetr::get_expressed_genes() per cell type, then union
    # This matches the NicheNet vignette approach:
    #   - Compute expressed genes for each sender cell type separately
    #   - Union across all sender cell types
    #   - A gene expressed in ANY sender cell type passes the threshold

    # Set identity to cell type column for get_expressed_genes()
    Seurat::Idents(object) <- object@meta.data[[group_col]]

    # Receiver expressed genes
    receiver_expressed <- nichenetr::get_expressed_genes(
        receiver, object, pct = min_pct, assay_oi = assay
    )

    # Sender expressed genes: per cell type, then union
    sender_expressed <- sender %>%
        unique() %>%
        lapply(nichenetr::get_expressed_genes, object,
               pct = min_pct, assay_oi = assay) %>%
        unlist() %>%
        unique()

    # Background: receiver expressed genes (matching NicheNet vignette)
    background_genes <- receiver_expressed

    list(
        sender_expressed = sender_expressed,
        receiver_expressed = receiver_expressed,
        background_genes = background_genes
    )
}


#' Get potential ligands from sender to receiver
#'
#' @description Identifies ligands expressed in sender cells that have
#'   receptors expressed in receiver cells.
#'
#' @param lr_network Ligand-receptor network data.frame
#' @param sender_expressed Genes expressed in sender
#' @param receiver_expressed Genes expressed in receiver
#'
#' @return Character vector of potential ligand names
#'
#' @keywords internal
.get_potential_ligands <- function(lr_network,
                                    sender_expressed,
                                    receiver_expressed) {
    # Determine column names (support both from/to and ligand/receptor)
    lig_col <- if ("ligand" %in% colnames(lr_network)) "ligand" else "from"
    rec_col <- if ("receptor" %in% colnames(lr_network)) "receptor" else "to"

    # Ligands expressed in sender
    ligands_sender <- unique(lr_network[[lig_col]][
        lr_network[[lig_col]] %in% sender_expressed])

    # Receptors expressed in receiver
    receptors_receiver <- unique(lr_network[[rec_col]][
        lr_network[[rec_col]] %in% receiver_expressed])

    # Ligands that have receptors in receiver
    lr_filtered <- lr_network[
        lr_network[[lig_col]] %in% ligands_sender &
        lr_network[[rec_col]] %in% receptors_receiver,
    ]

    unique(lr_filtered[[lig_col]])
}

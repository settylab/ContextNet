#' Export Functions
#'
#' @description
#' Functions for exporting ContextNet analysis results to various formats
#' for visualization with ccflowR (directly in R) and other tools like Cytoscape.
#'
#' @name export
NULL


#' Export network for Cytoscape visualization
#'
#' @description Export a ContextNetwork in GraphML or XGMML format for
#'   visualization in Cytoscape.
#'
#' @param network A ContextNetwork object
#' @param output_file Path to output file
#' @param format Export format: "graphml" (default) or "xgmml"
#' @param network_type Which network to export: "signaling" (default), "regulatory", or "both"
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return The output file path (invisibly)
#'
#' @examples
#' \dontrun{
#' network <- integrate_networks(grn)
#' export_to_cytoscape(network, "network.graphml")
#' }
#'
#' @export
export_to_cytoscape <- function(network,
                                 output_file,
                                 format = c("graphml", "xgmml"),
                                 network_type = c("signaling", "regulatory", "both"),
                                 verbose = TRUE) {
    # ── Validate inputs ──
    if (!inherits(network, "ContextNetwork")) {
        cli::cli_abort(c(
            "{.arg network} must be a ContextNetwork object",
            "i" = "Create one with {.fn integrate_networks}"
        ))
    }

    if (is.null(output_file) || !nzchar(output_file)) {
        cli::cli_abort("{.arg output_file} is required")
    }

    format <- match.arg(format)
    network_type <- match.arg(network_type)

    if (verbose) cli::cli_h2("Exporting to Cytoscape ({format})")

    # ── Get edges based on network type ──
    if (network_type == "signaling") {
        edges <- network@signaling_network
        graph_name <- "Signaling_Network"
    } else if (network_type == "regulatory") {
        edges <- network@regulatory_network
        graph_name <- "Regulatory_Network"
    } else {
        # Combine both networks
        sig <- network@signaling_network
        sig$network_type <- "signaling"

        reg <- network@regulatory_network
        reg$network_type <- "regulatory"

        edges <- rbind(sig, reg)
        graph_name <- "Combined_Network"
    }

    # Ensure required columns exist
    if (!"weight" %in% colnames(edges)) {
        edges$weight <- 1.0
    }

    if (verbose) {
        cli::cli_alert_info("Exporting {nrow(edges)} edges")
    }

    # ── Write output ──
    if (format == "graphml") {
        .write_graphml(edges, output_file, graph_name)
    } else {
        .write_xgmml(edges, output_file, graph_name)
    }

    if (verbose) {
        cli::cli_alert_success("Exported to {.path {output_file}}")
    }

    invisible(output_file)
}


#' Write GraphML format
#'
#' @param edges Edge data.frame with from, to, weight columns
#' @param output_file Output file path
#' @param graph_name Name for the graph
#'
#' @keywords internal
.write_graphml <- function(edges, output_file, graph_name) {
    # Get unique nodes
    nodes <- unique(c(edges$from, edges$to))

    # Build GraphML content
    lines <- c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns"',
        '         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"',
        '         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns',
        '         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">',
        '  <key id="weight" for="edge" attr.name="weight" attr.type="double"/>',
        '  <key id="label" for="node" attr.name="label" attr.type="string"/>'
    )

    # Add network_type key if present
    if ("network_type" %in% colnames(edges)) {
        lines <- c(lines,
            '  <key id="network_type" for="edge" attr.name="network_type" attr.type="string"/>'
        )
    }

    lines <- c(lines,
        paste0('  <graph id="', graph_name, '" edgedefault="directed">')
    )

    # Add nodes
    for (node in nodes) {
        lines <- c(lines,
            paste0('    <node id="', node, '">'),
            paste0('      <data key="label">', node, '</data>'),
            '    </node>'
        )
    }

    # Add edges
    for (i in seq_len(nrow(edges))) {
        edge_id <- paste0("e", i)
        lines <- c(lines,
            paste0('    <edge id="', edge_id, '" source="', edges$from[i], '" target="', edges$to[i], '">'),
            paste0('      <data key="weight">', round(edges$weight[i], 6), '</data>')
        )

        if ("network_type" %in% colnames(edges)) {
            lines <- c(lines,
                paste0('      <data key="network_type">', edges$network_type[i], '</data>')
            )
        }

        lines <- c(lines, '    </edge>')
    }

    lines <- c(lines,
        '  </graph>',
        '</graphml>'
    )

    writeLines(lines, output_file)
}


#' Write XGMML format
#'
#' @param edges Edge data.frame with from, to, weight columns
#' @param output_file Output file path
#' @param graph_name Name for the graph
#'
#' @keywords internal
.write_xgmml <- function(edges, output_file, graph_name) {
    # Get unique nodes
    nodes <- unique(c(edges$from, edges$to))

    # Build XGMML content
    lines <- c(
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graph xmlns="http://www.cs.rpi.edu/XGMML"',
        '       xmlns:dc="http://purl.org/dc/elements/1.1/"',
        '       xmlns:xlink="http://www.w3.org/1999/xlink"',
        '       xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"',
        '       xmlns:cy="http://www.cytoscape.org"',
        paste0('       label="', graph_name, '" directed="1">')
    )

    # Add nodes
    for (i in seq_along(nodes)) {
        node <- nodes[i]
        lines <- c(lines,
            paste0('  <node id="', i, '" label="', node, '">'),
            paste0('    <att name="name" type="string" value="', node, '"/>'),
            '  </node>'
        )
    }

    # Create node ID lookup
    node_ids <- setNames(seq_along(nodes), nodes)

    # Add edges
    for (i in seq_len(nrow(edges))) {
        source_id <- node_ids[edges$from[i]]
        target_id <- node_ids[edges$to[i]]
        weight <- round(edges$weight[i], 6)

        lines <- c(lines,
            paste0('  <edge source="', source_id, '" target="', target_id, '">'),
            paste0('    <att name="weight" type="real" value="', weight, '"/>'),
            paste0('    <att name="interaction" type="string" value="regulates"/>')
        )

        if ("network_type" %in% colnames(edges)) {
            lines <- c(lines,
                paste0('    <att name="network_type" type="string" value="', edges$network_type[i], '"/>')
            )
        }

        lines <- c(lines, '  </edge>')
    }

    lines <- c(lines, '</graph>')

    writeLines(lines, output_file)
}


#' Prepare data for ccflowR visualization
#'
#' @description Prepare ContextResult data in a format ready for direct use
#'   with the ccflowR R package for signaling flow visualization. ccflowR is
#'   an R package that can be used directly in R without needing to export
#'   files to disk.
#'
#' @param result A ContextResult object from analyze_communication()
#' @param expression_data Seurat object for computing z-score expression
#'   probabilities per cell type. Required for the full NicheNet weighted
#'   score pipeline. If NULL, falls back to normalized AUPR only.
#' @param network A ContextNetwork object for ligand-receptor interaction
#'   weights from the signaling network (lr_sig). If NULL, uses receptor_info
#'   from the result (which may lack weights).
#' @param sender Sender cell type(s) to include (default: from result parameters)
#' @param receiver Receiver cell type(s) to include (default: from result parameters)
#' @param top_n_ligands Maximum number of ligands to include (default: 50)
#' @param group_col Column name in Seurat metadata for cell type annotations
#'   (default: auto-detected from result parameters or common names)
#' @param excluded_celltypes Character vector of cell types to exclude from
#'   sender expansion (e.g., cell types with too few cells). Default: NULL.
#' @param min_cells Minimum number of cells required for a cell type to be
#'   included as a sender. Cell types with fewer than \code{min_cells} cells
#'   are automatically excluded. Requires \code{expression_data} to be a Seurat
#'   object. Set to NULL (default) to disable automatic exclusion. This is
#'   combined with \code{excluded_celltypes} (both are applied).
#' @param framework Output format: "nichenet" or "multinichenet" (default: "nichenet")
#' @param receptor_weight How to weight ligand-receptor edges for NicheNet results.
#'   \describe{
#'     \item{"network"}{(Default) Use L-R signaling network weight from NicheNet database.}
#'     \item{"expression"}{Use receptor gene expression (z-score probability in the
#'       receiver cell type). Requires \code{expression_data}. This mirrors the
#'       MultiNicheNet approach where receptor expression drives edge weighting.}
#'   }
#'
#' @return A data.frame with per-receiver raw scores (7 columns):
#'   \describe{
#'     \item{sender}{Sender cell type (one row per sender cell type)}
#'     \item{ligand}{Ligand gene name}
#'     \item{receptor}{Receptor gene name}
#'     \item{receiver}{Receiver cell type}
#'     \item{weighted_score}{(zscore_prob + normalized_aupr) / 2 per sender}
#'     \item{ligand_expression}{Raw z-score expression probability (before combining with AUPR)}
#'     \item{weight}{L-R interaction weight from signaling network}
#'   }
#'   To get ccflowR-ready output, pass one or more of these data.frames to
#'   \code{\link{combine_sankey_data}()} which adds global scaling and thresholding.
#'
#' @details
#' This is the recommended way to visualize ContextNet results with ccflowR.
#' The ccflowR package is an R package that provides sankey/river plots,
#' heatmaps, and other signaling flow visualizations directly in R.
#'
#' ## NicheNet Framework Pipeline
#'
#' The pipeline replicates the manual NicheNet Sankey workflow:
#' 1. Z-score each gene across ALL cells, then average per cell type
#' 2. Convert averaged z-scores to probability via pnorm()
#' 3. Normalize AUPR to 0-1 (positive AUPR only, min-max)
#' 4. Combine: (zscore_prob + normalized_aupr) / 2 per sender cell type
#' 5. Join with L-R interaction weights from signaling network
#'
#' The output contains raw scores only. Use \code{\link{combine_sankey_data}()}
#' to combine multiple receivers, apply global min-max scaling, and threshold.
#'
#' ## MultiNicheNet Framework Pipeline
#'
#' For MultiNicheNet results, the function extracts the
#' \code{group_prioritization_tbl} stored in the ContextResult, filters to the
#' condition of interest, and maps prioritization columns to the same 7-column
#' format. Since MultiNicheNet analyzes all receivers in one run, call
#' \code{prepare_for_ccflowr()} once and pipe directly to
#' \code{\link{combine_sankey_data}()} for scaling and thresholding.
#'
#' @examples
#' \dontrun{
#' # --- NicheNet workflow (one call per receiver) ---
#' data_basal <- prepare_for_ccflowr(result_basal, seurat_obj, network)
#' data_parabasal <- prepare_for_ccflowr(result_parabasal, seurat_obj, network)
#'
#' # Combine, scale, and threshold in one step
#' sankey_data <- combine_sankey_data(data_basal, data_parabasal,
#'                                     threshold_quantile = 0.95)
#'
#' # Use directly with ccflowR in R (no file export needed)
#' library(ccflowR)
#' sankey_flow(sankey_data, cell_colors)
#'
#' # --- MultiNicheNet workflow (single call, all receivers) ---
#' mn_data <- prepare_for_ccflowr(result_multi, framework = "multinichenet")
#' sankey_data <- combine_sankey_data(mn_data, threshold_quantile = 0.95)
#' sankey_flow(sankey_data, cell_colors)
#' }
#'
#' @export
prepare_for_ccflowr <- function(result,
                                 expression_data = NULL,
                                 network = NULL,
                                 sender = NULL,
                                 receiver = NULL,
                                 top_n_ligands = 50,
                                 group_col = NULL,
                                 excluded_celltypes = NULL,
                                 min_cells = NULL,
                                 framework = c("nichenet", "multinichenet"),
                                 receptor_weight = c("network", "expression")) {
    # ── Validate inputs ──
    if (!inherits(result, "ContextResult")) {
        cli::cli_abort(c(
            "{.arg result} must be a ContextResult object",
            "i" = "Create one with {.fn analyze_communication}"
        ))
    }

    # Auto-detect framework from result when not explicitly specified
    if (missing(framework) && !is.null(result@analysis_type) &&
        result@analysis_type == "multi_condition") {
        framework <- "multinichenet"
    } else {
        framework <- match.arg(framework)
    }

    # Use parameters from result if not specified
    if (is.null(sender)) {
        sender <- result@parameters$sender
    }
    if (is.null(receiver)) {
        receiver <- result@parameters$receiver
    }
    if (is.null(group_col)) {
        group_col <- result@parameters$group_col %||%
            result@parameters$ident_col
    }

    # ── Validate min_cells requires expression_data ──
    if (!is.null(min_cells) && min_cells > 0 && is.null(expression_data)) {
        cli::cli_warn(c(
            "{.arg min_cells} requires {.arg expression_data} (a Seurat object) to count cells per cell type",
            "i" = "Ignoring {.arg min_cells} because no expression data was provided"
        ))
        min_cells <- NULL
    }

    # ── Check receiver min_cells ──
    if (!is.null(min_cells) && min_cells > 0 &&
        !is.null(expression_data) && inherits(expression_data, "Seurat") &&
        !is.null(group_col) && !is.null(receiver)) {
        ct_counts <- table(expression_data@meta.data[[group_col]])
        for (rcv in receiver) {
            if (rcv %in% names(ct_counts) && ct_counts[[rcv]] < min_cells) {
                cli::cli_warn(c(
                    "Receiver {.val {rcv}} has only {ct_counts[[rcv]]} cells (< {.arg min_cells}={min_cells})",
                    "i" = "Results may be unreliable with too few receiver cells"
                ))
            }
        }
    }

    # ── Get top ligands ──
    top_ligs <- utils::head(result@top_ligands, top_n_ligands)

    # ── Get ligand activities ──
    all_la <- result@ligand_activities
    la <- all_la

    # Handle different column names
    ligand_col <- if ("test_ligand" %in% colnames(la)) {
        "test_ligand"
    } else if ("ligand" %in% colnames(la)) {
        "ligand"
    } else {
        NULL
    }

    if (!is.null(ligand_col) && nrow(la) > 0) {
        la <- la[la[[ligand_col]] %in% top_ligs, , drop = FALSE]
    }

    receptor_weight <- match.arg(receptor_weight)

    # ── Framework-specific processing ──
    if (framework == "nichenet") {
        # Pass all_la so AUPR normalization uses the full range
        # (matching the manual NicheNet notebook workflow)
        out <- .prepare_nichenet_sankey(result, la, top_ligs, sender, receiver,
                                        expression_data, network, group_col,
                                        excluded_celltypes, all_la, min_cells,
                                        receptor_weight)
    } else {
        out <- .prepare_multinichenet_sankey(result, la, top_ligs, sender, receiver)
    }

    # Replace NAs with small positive value
    numeric_cols <- c("weighted_score", "weight")
    for (col in intersect(numeric_cols, colnames(out))) {
        out[[col]][is.na(out[[col]])] <- 0.01
    }

    # Remove duplicates
    out <- unique(out)

    # Set framework attribute
    attr(out, "framework") <- framework

    cli::cli_alert_success("Prepared {nrow(out)} links for ccflowR visualization")

    out
}


#' Combine and scale sankey data from multiple receivers
#'
#' @description Combines per-receiver data.frames from
#'   \code{\link{prepare_for_ccflowr}()}, applies global min-max scaling, maps
#'   ccflowR-compatible columns, and optionally filters by a quantile threshold.
#'
#' @param ... One or more data.frames returned by \code{prepare_for_ccflowr()},
#'   or a single list of such data.frames.
#' @param threshold_quantile Quantile threshold for filtering rows (default:
#'   0.95). Set to 0 to skip thresholding and return all rows with scaled
#'   columns.
#'
#' @return A data.frame with 11 columns: the 7 raw columns from
#'   \code{prepare_for_ccflowr()} plus 4 scaled columns:
#'   \describe{
#'     \item{sender}{Sender cell type}
#'     \item{ligand}{Ligand gene name}
#'     \item{receptor}{Receptor gene name}
#'     \item{receiver}{Receiver cell type}
#'     \item{weighted_score}{Raw combined score per sender}
#'     \item{ligand_expression}{Raw z-score expression probability}
#'     \item{weight}{L-R interaction weight from signaling network}
#'     \item{scaled_weighted_score}{Global min-max scaled weighted_score}
#'     \item{scaled_activity}{Same as scaled_weighted_score (ccflowR compat)}
#'     \item{scaled_ligand_expression}{Min-max scaled ligand_expression (ccflowR compat)}
#'     \item{scaled_receptor_expression}{Min-max scaled L-R weight (ccflowR compat)}
#'   }
#'
#' @details
#' This function encapsulates the combine-rescale-threshold pattern that is
#' needed when analyzing multiple receiver cell types. Each call to
#' \code{prepare_for_ccflowr()} produces per-receiver raw scores;
#' \code{combine_sankey_data()} then:
#' \enumerate{
#'   \item Combines all inputs via \code{rbind()}
#'   \item Applies global min-max scaling to \code{weighted_score}
#'   \item Maps ccflowR-compatible column aliases
#'   \item Filters to rows at or above the \code{threshold_quantile}
#' }
#'
#' Works with both NicheNet and MultiNicheNet framework results from
#' \code{prepare_for_ccflowr()}.
#'
#' @examples
#' \dontrun{
#' # Prepare per-receiver data
#' data_basal <- prepare_for_ccflowr(result_basal, seurat_obj, network)
#' data_parabasal <- prepare_for_ccflowr(result_parabasal, seurat_obj, network)
#' data_suprabasal <- prepare_for_ccflowr(result_suprabasal, seurat_obj, network)
#'
#' # Combine, scale globally, and threshold at 95th percentile
#' sankey_data <- combine_sankey_data(
#'     data_basal, data_parabasal, data_suprabasal,
#'     threshold_quantile = 0.95
#' )
#'
#' # Or get all data without thresholding
#' full_data <- combine_sankey_data(
#'     data_basal, data_parabasal, data_suprabasal,
#'     threshold_quantile = 0
#' )
#'
#' # Also accepts a list
#' sankey_data <- combine_sankey_data(
#'     list(data_basal, data_parabasal),
#'     threshold_quantile = 0.95
#' )
#' }
#'
#' @export
combine_sankey_data <- function(..., threshold_quantile = 0.95) {
    inputs <- list(...)

    # Accept a single list of data.frames
    if (length(inputs) == 1 && is.list(inputs[[1]]) &&
        !is.data.frame(inputs[[1]])) {
        inputs <- inputs[[1]]
    }

    if (length(inputs) == 0) {
        cli::cli_abort("At least one data.frame is required")
    }

    # Validate required columns
    required_cols <- c("sender", "ligand", "receptor", "receiver",
                       "weighted_score", "weight")
    for (i in seq_along(inputs)) {
        if (!is.data.frame(inputs[[i]])) {
            cli::cli_abort("Input {i} is not a data.frame")
        }
        missing <- setdiff(required_cols, colnames(inputs[[i]]))
        if (length(missing) > 0) {
            cli::cli_abort(c(
                "Input {i} is missing required columns: {.val {missing}}",
                "i" = "Use {.fn prepare_for_ccflowr} to generate the input data"
            ))
        }
    }

    # Strip any pre-existing scaled_* columns (handles legacy data)
    scaled_cols <- c("scaled_weighted_score", "scaled_activity",
                     "scaled_ligand_expression", "scaled_receptor_expression")
    inputs <- lapply(inputs, function(df) {
        drop <- intersect(scaled_cols, colnames(df))
        if (length(drop) > 0) {
            df[, setdiff(colnames(df), drop), drop = FALSE]
        } else {
            df
        }
    })

    # Ensure ligand_expression column exists for rbind compatibility
    inputs <- lapply(inputs, function(df) {
        if (!("ligand_expression" %in% colnames(df))) {
            df$ligand_expression <- df$weighted_score
        }
        df
    })

    # Combine
    combined <- do.call(rbind, inputs)

    if (nrow(combined) == 0) {
        # Return empty with all 10 columns
        empty <- .empty_sankey_df()
        empty$scaled_weighted_score <- numeric()
        empty$scaled_activity <- numeric()
        empty$scaled_ligand_expression <- numeric()
        empty$scaled_receptor_expression <- numeric()
        return(empty)
    }

    # Global min-max scale weighted_score
    combined$scaled_weighted_score <- .scale_01(combined$weighted_score)

    # Map ccflowR-compatible columns
    combined$scaled_activity <- combined$scaled_weighted_score
    # Use pure ligand expression (zscore_prob) for sender->ligand edge weighting
    combined$scaled_ligand_expression <- pmax(.scale_01(combined$ligand_expression), 1e-4)
    combined$scaled_receptor_expression <- pmax(.scale_01(combined$weight), 1e-4)

    # Apply threshold
    if (threshold_quantile > 0) {
        threshold <- stats::quantile(combined$scaled_weighted_score,
                                     threshold_quantile, na.rm = TRUE)
        combined <- combined[combined$scaled_weighted_score >= threshold, ,
                             drop = FALSE]
    }

    # Reorder columns
    col_order <- c("sender", "ligand", "receptor", "receiver",
                   "weighted_score", "ligand_expression", "weight",
                   "scaled_weighted_score", "scaled_activity",
                   "scaled_ligand_expression", "scaled_receptor_expression")
    combined <- combined[, intersect(col_order, colnames(combined)), drop = FALSE]

    combined
}


#' Prepare NicheNet-style sankey data
#'
#' @description Replicates the manual NicheNet Sankey pipeline:
#'   1. Z-score all cells, average per cell type, pnorm
#'   2. Normalize AUPR (positive only, min-max)
#'   3. Combine: (zscore_prob + normalized_aupr) / 2 per sender
#'   4. Join with L-R weights from signaling network
#'
#' @param result ContextResult object
#' @param la Filtered ligand activities data.frame (top N ligands)
#' @param top_ligs Top ligands to include
#' @param sender Sender cell type(s)
#' @param receiver Receiver cell type(s)
#' @param expression_data Seurat object (or NULL)
#' @param network ContextNetwork object (or NULL)
#' @param group_col Cell type column name
#' @param excluded_celltypes Cell types to exclude from sender expansion
#' @param all_la All ligand activities (for AUPR normalization range)
#' @param min_cells Minimum cell count threshold for sender cell types (or NULL)
#'
#' @return data.frame with sankey columns
#' @keywords internal
.prepare_nichenet_sankey <- function(result, la, top_ligs, sender, receiver,
                                     expression_data, network, group_col,
                                     excluded_celltypes, all_la = NULL,
                                     min_cells = NULL,
                                     receptor_weight = "network") {
    if (nrow(la) == 0) {
        return(.empty_sankey_df())
    }

    # ── Standardize ligand column ──
    ligand_col <- if ("test_ligand" %in% colnames(la)) "test_ligand" else "ligand"
    if (ligand_col != "ligand" && !("ligand" %in% colnames(la))) {
        la$ligand <- la[[ligand_col]]
    }

    # ── Step 1: Normalize AUPR (positive only, min-max) ──
    # Use ALL ligand activities for the normalization range (not just top N).
    # This matches the manual NicheNet notebook workflow where normalized_score
    # is computed across all positive-AUPR ligands, giving top ligands scores
    # closer to 1.
    if (nrow(la) > 0 && "aupr_corrected" %in% colnames(la)) {
        # Determine normalization range from all ligands (if available)
        norm_source <- if (!is.null(all_la) && nrow(all_la) > 0 &&
                           "aupr_corrected" %in% colnames(all_la)) {
            all_la
        } else {
            la
        }
        norm_pos <- norm_source[norm_source$aupr_corrected > 0, , drop = FALSE]

        if (nrow(norm_pos) > 0) {
            max_aupr <- max(norm_pos$aupr_corrected, na.rm = TRUE)
            min_aupr <- min(norm_pos$aupr_corrected, na.rm = TRUE)
        } else {
            max_aupr <- min_aupr <- 0
        }

        # Apply normalization to the top-N ligands only
        la_pos <- la[la$aupr_corrected > 0, , drop = FALSE]
        if (nrow(la_pos) > 0 && max_aupr > min_aupr) {
            la_pos$normalized_score <- (la_pos$aupr_corrected - min_aupr) /
                (max_aupr - min_aupr)
            la <- la_pos
        } else if (nrow(la_pos) > 0) {
            la_pos$normalized_score <- 1
            la <- la_pos
        } else {
            la$normalized_score <- 1
        }
    } else {
        la$normalized_score <- 1
    }

    # ── Step 2: Compute z-score expression probability per cell type ──
    zscore_prob <- NULL
    if (!is.null(expression_data)) {
        zscore_prob <- .compute_zscore_prob(expression_data, la$ligand, group_col)
    }

    if (is.null(zscore_prob)) {
        # No expression data: fall back to normalized AUPR only
        return(.prepare_nichenet_sankey_fallback(la, sender, receiver,
                                                 result, network, top_ligs))
    }

    # ── Step 3: Combine zscore_prob with normalized AUPR ──
    # For each ligand, for each cell type: (zscore_prob + normalized_score) / 2
    ct_cols <- setdiff(colnames(zscore_prob), "Gene")
    df_ws <- zscore_prob[zscore_prob$Gene %in% la$ligand, , drop = FALSE]
    df_ws <- merge(df_ws,
                   la[, c("ligand", "normalized_score"), drop = FALSE],
                   by.x = "Gene", by.y = "ligand", all.x = TRUE)

    # Save raw expression probabilities before combining with activity
    raw_expr <- df_ws[, ct_cols, drop = FALSE]

    for (ct in ct_cols) {
        df_ws[[ct]] <- (df_ws[[ct]] + df_ws$normalized_score) / 2
    }
    df_ws$normalized_score <- NULL
    colnames(df_ws)[colnames(df_ws) == "Gene"] <- "ligand"

    # ── Step 4: Pivot to long (one row per ligand x celltype) ──
    rows <- lapply(ct_cols, function(ct) {
        data.frame(
            ligand = df_ws$ligand,
            celltype = ct,
            weighted_score = df_ws[[ct]],
            ligand_expression = raw_expr[[ct]],
            stringsAsFactors = FALSE
        )
    })
    df_long <- do.call(rbind, rows)

    # ── Step 5: Exclude cell types ──
    # Do NOT filter to sender cell types here. The manual NicheNet notebook
    # includes ALL cell types as potential senders (only excluding specific ones
    # like "dividing" or those with < 500 cells). Filtering to sender would
    # remove minor cell types (e.g., B Cell, Melanocyte) that contribute to the
    # score distribution and affect the 95th percentile threshold.

    # Auto-exclude cell types with fewer than min_cells cells
    if (!is.null(min_cells) && min_cells > 0 &&
        !is.null(expression_data) && inherits(expression_data, "Seurat") &&
        !is.null(group_col)) {
        ct_counts <- table(expression_data@meta.data[[group_col]])
        low_ct <- names(ct_counts)[ct_counts < min_cells]
        if (length(low_ct) > 0) {
            ct_info <- paste0(low_ct, " (", ct_counts[low_ct], ")")
            cli::cli_alert_info(
                "Excluding {length(low_ct)} cell type{?s} with <{min_cells} cells: {.val {ct_info}}"
            )
        }
        excluded_celltypes <- unique(c(excluded_celltypes, low_ct))
    }

    if (!is.null(excluded_celltypes) && length(excluded_celltypes) > 0) {
        n_before <- nrow(df_long)
        df_long <- df_long[!df_long$celltype %in% excluded_celltypes, , drop = FALSE]
        cli::cli_alert_info(
            "Kept {length(unique(df_long$celltype))} sender cell types ({nrow(df_long)}/{n_before} rows)"
        )
    }

    # ── Step 6: Get L-R pairs with weights ──
    lr_df <- .get_lr_weights(result, network, la$ligand, top_ligs)

    if (nrow(lr_df) == 0) {
        return(.empty_sankey_df())
    }

    # ── Step 6b: Replace network weights with receptor expression if requested ──
    if (receptor_weight == "expression" && !is.null(expression_data)) {
        receptor_genes <- unique(lr_df$receptor)
        receptor_zscore <- .compute_zscore_prob(expression_data, receptor_genes,
                                                group_col)
        if (!is.null(receptor_zscore) && receiver[1] %in% colnames(receptor_zscore)) {
            # Look up each receptor's z-score probability in the receiver cell type
            receptor_expr_map <- stats::setNames(
                receptor_zscore[[receiver[1]]],
                receptor_zscore$Gene
            )
            lr_df$weight <- vapply(lr_df$receptor, function(r) {
                if (r %in% names(receptor_expr_map)) receptor_expr_map[[r]] else 0.5
            }, numeric(1))
            cli::cli_alert_info(
                "Using receptor expression (z-score probability in {.val {receiver[1]}}) for L-R weighting"
            )
        } else {
            cli::cli_warn(
                "Cannot compute receptor expression; falling back to network weights"
            )
        }
    }

    # ── Step 7: Join L-R pairs (many-to-many: ligand x celltype x receptor) ──
    full_data <- merge(df_long, lr_df, by = "ligand")
    full_data$receiver <- receiver[1]

    # ── Step 8: Rename celltype -> sender ──
    colnames(full_data)[colnames(full_data) == "celltype"] <- "sender"

    # ── Reorder columns (raw scores only; use combine_sankey_data() for scaling) ──
    col_order <- c("sender", "ligand", "receptor", "receiver",
                   "weighted_score", "ligand_expression", "weight")
    full_data <- full_data[, intersect(col_order, colnames(full_data)), drop = FALSE]

    full_data
}


#' Get L-R pairs with weights from signaling network
#'
#' @param result ContextResult object
#' @param network ContextNetwork object (or NULL)
#' @param ligands Ligands to include
#' @param top_ligs Top ligands from analysis
#'
#' @return data.frame with ligand, receptor, weight columns
#' @keywords internal
.get_lr_weights <- function(result, network, ligands, top_ligs) {
    # Get receptor_info pairs from result
    ri <- result@receptor_info
    ri_pairs <- data.frame(ligand = character(), receptor = character(),
                           stringsAsFactors = FALSE)

    if (nrow(ri) > 0) {
        if ("from" %in% colnames(ri) && "to" %in% colnames(ri)) {
            keep_cols <- intersect(c("from", "to", "weight"), colnames(ri))
            ri_pairs <- ri[, keep_cols, drop = FALSE]
            colnames(ri_pairs)[colnames(ri_pairs) == "from"] <- "ligand"
            colnames(ri_pairs)[colnames(ri_pairs) == "to"] <- "receptor"
        } else if ("ligand" %in% colnames(ri) && "receptor" %in% colnames(ri)) {
            keep_cols <- intersect(c("ligand", "receptor", "weight"), colnames(ri))
            ri_pairs <- ri[, keep_cols, drop = FALSE]
        }
        ri_pairs <- ri_pairs[ri_pairs$ligand %in% top_ligs, , drop = FALSE]
    }

    # Deduplicate L-R pairs (lr_network has multiple sources per pair)
    ri_pairs <- unique(ri_pairs[, c("ligand", "receptor"), drop = FALSE])

    # Try to get weights from ContextNetwork signaling_network (= lr_sig)
    if (!is.null(network) && inherits(network, "ContextNetwork") &&
        nrow(network@signaling_network) > 0) {
        lr_sig <- network@signaling_network
        lr_sig_w <- data.frame(
            ligand = lr_sig$from,
            receptor = lr_sig$to,
            weight = lr_sig$weight,
            stringsAsFactors = FALSE
        )
        lr_sig_w <- lr_sig_w[!duplicated(lr_sig_w[, c("ligand", "receptor")]), ,
                             drop = FALSE]

        if (nrow(ri_pairs) > 0) {
            # Inner join: only keep L-R pairs present in signaling network
            # (matches nichenetr::get_weighted_ligand_receptor_links behavior)
            lr_df <- merge(ri_pairs, lr_sig_w,
                           by = c("ligand", "receptor"))
        } else {
            # No receptor_info; filter lr_sig for top ligands
            lr_df <- lr_sig_w[lr_sig_w$ligand %in% top_ligs, , drop = FALSE]
        }
        return(lr_df)
    }

    # No network: use receptor_info with existing weight or default
    if (nrow(ri_pairs) > 0) {
        lr_df <- ri_pairs
        if (!"weight" %in% colnames(lr_df)) {
            lr_df$weight <- 1
        }
        return(lr_df)
    }

    # Fallback: use ligand_target_links as pseudo L-R data
    ltl <- result@ligand_target_links
    if (nrow(ltl) > 0) {
        ltl <- ltl[ltl$ligand %in% top_ligs, , drop = FALSE]
        lr_df <- data.frame(
            ligand = ltl$ligand,
            receptor = ltl$target,
            weight = if ("weight" %in% colnames(ltl)) ltl$weight else
                     if ("score" %in% colnames(ltl)) ltl$score else 1,
            stringsAsFactors = FALSE
        )
        return(lr_df)
    }

    data.frame(ligand = character(), receptor = character(),
               weight = numeric(), stringsAsFactors = FALSE)
}


#' Fallback NicheNet sankey when no expression data is available
#'
#' @description Uses normalized AUPR only (no z-score), with sender[1].
#' @keywords internal
.prepare_nichenet_sankey_fallback <- function(la, sender, receiver,
                                              result, network, top_ligs) {
    lr_df <- .get_lr_weights(result, network, la$ligand, top_ligs)

    if (nrow(lr_df) == 0 || nrow(la) == 0) {
        return(.empty_sankey_df())
    }

    # Use normalized_score directly as weighted_score
    la_slim <- la[, c("ligand", "normalized_score"), drop = FALSE]
    colnames(la_slim)[colnames(la_slim) == "normalized_score"] <- "weighted_score"

    out <- merge(la_slim, lr_df, by = "ligand")
    out$sender <- sender[1]
    out$receiver <- receiver[1]
    out$ligand_expression <- out$weighted_score  # no raw expression available

    col_order <- c("sender", "ligand", "receptor", "receiver",
                   "weighted_score", "ligand_expression", "weight")
    out <- out[, intersect(col_order, colnames(out)), drop = FALSE]

    out
}


#' Empty sankey data.frame with correct structure
#' @keywords internal
.empty_sankey_df <- function() {
    data.frame(
        sender = character(), ligand = character(),
        receptor = character(), receiver = character(),
        weighted_score = numeric(), ligand_expression = numeric(),
        weight = numeric(),
        stringsAsFactors = FALSE
    )
}


#' Prepare MultiNicheNet-style sankey data
#'
#' @description Maps MultiNicheNet output to sankey columns:
#'   - scaled_activity: prioritization_score
#'   - scaled_ligand_expression: scaled_lfc_ligand
#'   - scaled_receptor_expression: scaled_lfc_receptor
#'
#' @param result ContextResult object
#' @param la Filtered ligand activities data.frame
#' @param top_ligs Top ligands to include
#' @param sender Sender cell type(s)
#' @param receiver Receiver cell type(s)
#'
#' @return data.frame with sankey columns
#' @keywords internal
.prepare_multinichenet_sankey <- function(result, la, top_ligs, sender, receiver) {

    # ── Try to use group_prioritization_tbl (stored by predict_ligand_activity_multinichenet) ──
    gpt <- result@expression_data$group_prioritization_tbl

    if (!is.null(gpt) && is.data.frame(gpt) && nrow(gpt) > 0 &&
        all(c("sender", "receiver", "ligand", "receptor",
              "prioritization_score", "group", "contrast") %in% colnames(gpt))) {

        # Filter to condition of interest and the forward contrast
        condition_oi <- result@parameters$condition_oi %||%
            result@parameters$ident_1
        condition_ref <- result@parameters$ident_2
        contrast_oi <- paste0(condition_oi, "-", condition_ref)

        prio <- gpt[gpt$group == condition_oi & gpt$contrast == contrast_oi, ,
                     drop = FALSE]

        if (nrow(prio) == 0) {
            # Try without contrast filter (older multinichenetr versions)
            prio <- gpt[gpt$group == condition_oi, , drop = FALSE]
        }

        # Filter to top ligands
        prio <- prio[prio$ligand %in% top_ligs, , drop = FALSE]

        if (nrow(prio) > 0) {
            # Map to the same 7-column format as .prepare_nichenet_sankey():
            #   sender, ligand, receptor, receiver, weighted_score, ligand_expression, weight
            has_lfc_lig <- "scaled_lfc_ligand" %in% colnames(prio)
            has_lfc_rec <- "scaled_lfc_receptor" %in% colnames(prio)

            out <- data.frame(
                sender = prio$sender,
                ligand = prio$ligand,
                receptor = prio$receptor,
                receiver = prio$receiver,
                weighted_score = prio$prioritization_score,
                ligand_expression = if (has_lfc_lig) prio$scaled_lfc_ligand
                                    else prio$prioritization_score,
                weight = if (has_lfc_rec) prio$scaled_lfc_receptor
                         else prio$prioritization_score,
                stringsAsFactors = FALSE
            )
            return(out)
        }
    }

    # ── Fallback: no group_prioritization_tbl available ──
    # Use ligand_activities if they have the needed columns
    if (nrow(la) > 0 && "aupr_corrected" %in% colnames(la)) {
        # Get receptor info for L-R pairing
        ri <- result@receptor_info
        if (nrow(ri) > 0 && all(c("from", "to") %in% colnames(ri))) {
            ri_pairs <- unique(ri[ri$from %in% top_ligs,
                                  c("from", "to"), drop = FALSE])
            colnames(ri_pairs) <- c("ligand", "receptor")
        } else {
            ri_pairs <- data.frame(ligand = top_ligs, receptor = top_ligs,
                                   stringsAsFactors = FALSE)
        }

        # Merge activities with L-R pairs
        la_sub <- la[la$test_ligand %in% top_ligs, , drop = FALSE]
        la_sub <- la_sub[!duplicated(la_sub$test_ligand), , drop = FALSE]
        la_scores <- data.frame(
            ligand = la_sub$test_ligand,
            score = la_sub$aupr_corrected,
            stringsAsFactors = FALSE
        )

        out <- merge(ri_pairs, la_scores, by = "ligand", all.x = TRUE)
        out$score[is.na(out$score)] <- 0

        # Expand to all senders
        rows <- lapply(sender, function(s) {
            data.frame(
                sender = s,
                ligand = out$ligand,
                receptor = out$receptor,
                receiver = receiver[1],
                weighted_score = out$score,
                ligand_expression = out$score,
                weight = out$score,
                stringsAsFactors = FALSE
            )
        })
        return(do.call(rbind, rows))
    }

    .empty_sankey_df()
}


#' Compute z-score expression probability per cell type
#'
#' @description Replicates the manual NicheNet z-score pipeline:
#'   1. Z-score each gene across ALL cells (equivalent to ScaleData)
#'   2. Average z-score per gene per cell type
#'   3. Convert to probability via pnorm()
#'
#' @param expression_data Seurat object
#' @param genes Genes to compute probability for
#' @param group_col Column name in Seurat metadata for cell type annotations
#'
#' @return data.frame with Gene column and one column per cell type (pnorm values),
#'   or NULL if computation fails
#' @keywords internal
.compute_zscore_prob <- function(expression_data, genes, group_col = NULL) {
    if (!inherits(expression_data, "Seurat")) {
        return(NULL)
    }
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        return(NULL)
    }

    # ── Detect group column ──
    if (is.null(group_col)) {
        for (col in c("cell_type", "celltype", "CellType", "Fine_Cell_Types",
                       "seurat_clusters")) {
            if (col %in% colnames(expression_data@meta.data)) {
                group_col <- col
                break
            }
        }
    }
    if (is.null(group_col) ||
        !group_col %in% colnames(expression_data@meta.data)) {
        cli::cli_warn("Cannot determine cell type column from Seurat metadata")
        return(NULL)
    }

    # ── Get normalized expression (data layer) ──
    expr_mat <- Seurat::GetAssayData(expression_data, layer = "data")
    genes_avail <- intersect(genes, rownames(expr_mat))
    if (length(genes_avail) == 0) {
        return(NULL)
    }

    # Convert to dense matrix for subsetting
    expr_sub <- as.matrix(expr_mat[genes_avail, , drop = FALSE])  # genes x cells

    # ── Z-score each gene across ALL cells (same as ScaleData) ──
    gene_means <- rowMeans(expr_sub)
    gene_sds <- apply(expr_sub, 1, stats::sd)
    gene_sds[gene_sds == 0] <- 1  # avoid division by zero
    z_mat <- (expr_sub - gene_means) / gene_sds
    # Clip to [-10, 10] to match Seurat::ScaleData default (scale.max = 10)
    z_mat[z_mat > 10] <- 10
    z_mat[z_mat < -10] <- -10

    # ── Average z-score per gene per cell type ──
    cell_types <- as.character(expression_data@meta.data[[group_col]])
    unique_cts <- unique(cell_types)

    avg_z <- sapply(unique_cts, function(ct) {
        cells <- which(cell_types == ct)
        if (length(cells) > 0) {
            rowMeans(z_mat[, cells, drop = FALSE], na.rm = TRUE)
        } else {
            rep(0, nrow(z_mat))
        }
    })

    # Ensure avg_z is a matrix (sapply may return vector for single cell type)
    if (!is.matrix(avg_z)) {
        avg_z <- matrix(avg_z, nrow = length(genes_avail),
                        dimnames = list(genes_avail, unique_cts))
    }

    # ── Convert to probability via pnorm ──
    avg_z_prob <- stats::pnorm(avg_z)

    # ── Return as data.frame: Gene x cell_types ──
    result <- data.frame(Gene = rownames(avg_z_prob), check.names = FALSE,
                         stringsAsFactors = FALSE)
    for (ct in unique_cts) {
        result[[ct]] <- avg_z_prob[, ct]
    }

    result
}


#' Scale values to 0-1 range (min-max normalization)
#' @param x Numeric vector
#' @return Scaled numeric vector
#' @keywords internal
.scale_01 <- function(x) {
    x <- as.numeric(x)
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0.5, length(x)))
    (x - rng[1]) / diff(rng)
}


#' Extract expression data from Seurat or matrix
#' @param data Seurat object or expression matrix
#' @param sender Sender cell types
#' @param receiver Receiver cell types
#' @return List with ligand_expr and receptor_expr
#' @keywords internal
.extract_expression <- function(data, sender, receiver) {
    result <- list(ligand_expr = NULL, receptor_expr = NULL)

    if (inherits(data, "Seurat")) {
        if (!requireNamespace("Seurat", quietly = TRUE)) {
            return(result)
        }

        # Get expression matrix (use layer for Seurat v5)
        expr_mat <- Seurat::GetAssayData(data, layer = "data")

        # Get cell type column
        if ("cell_type" %in% colnames(data@meta.data)) {
            ct_col <- "cell_type"
        } else if ("celltype" %in% colnames(data@meta.data)) {
            ct_col <- "celltype"
        } else {
            return(result)
        }

        # Get sender cells
        sender_cells <- colnames(data)[data@meta.data[[ct_col]] %in% sender]
        if (length(sender_cells) > 0) {
            result$ligand_expr <- Matrix::rowMeans(expr_mat[, sender_cells, drop = FALSE])
        }

        # Get receiver cells
        receiver_cells <- colnames(data)[data@meta.data[[ct_col]] %in% receiver]
        if (length(receiver_cells) > 0) {
            result$receptor_expr <- Matrix::rowMeans(expr_mat[, receiver_cells, drop = FALSE])
        }
    } else if (is.matrix(data) || inherits(data, "Matrix")) {
        # Assume rows are genes, use row means
        result$ligand_expr <- rowMeans(data)
        result$receptor_expr <- rowMeans(data)
    }

    result
}


#' Export ContextResult to RDS file
#'
#' @description Simple export of a ContextResult object as an RDS file
#'   for later reloading in R.
#'
#' @param result A ContextResult object
#' @param output_file Path to output RDS file
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return The output file path (invisibly)
#'
#' @examples
#' \dontrun{
#' result <- analyze_communication(seurat_obj, network, ...)
#' export_result(result, "analysis_result.rds")
#' }
#'
#' @export
export_result <- function(result, output_file, verbose = TRUE) {
    # ── Validate inputs ──
    if (!inherits(result, "ContextResult")) {
        cli::cli_abort(c(
            "{.arg result} must be a ContextResult object",
            "i" = "Create one with {.fn analyze_communication}"
        ))
    }

    if (is.null(output_file) || !nzchar(output_file)) {
        cli::cli_abort("{.arg output_file} is required")
    }

    # ── Save ──
    saveRDS(result, output_file)

    if (verbose) {
        cli::cli_alert_success("Exported ContextResult to {.path {output_file}}")
    }

    invisible(output_file)
}

#' Parameter Tuning for Network Integration
#'
#' @description
#' Functions for tuning damping_factor and ltf_cutoff parameters when
#' constructing ligand-target matrices with context-specific GRNs.
#'
#' @name tune-params
NULL


#' Tune integration parameters based on GRN propagation consistency
#'
#' @description
#' Find optimal damping_factor and ltf_cutoff by maximizing how well
#' the ligand-target matrix preserves TF-target relationships from the GRN.
#'
#' The logic: If a TF is downstream of a ligand (via the signaling network),
#' then the TF's targets (from your GRN) should be predicted as targets of
#' that ligand in the final ligand-target matrix.
#'
#' @param grn A ContextGRN object or data.frame with columns: from (TF), to (target), weight
#' @param lr_network Ligand-receptor network (data.frame with from, to columns)
#' @param weighted_networks List with lr_sig and gr components (from NicheNet)
#' @param damping_range Numeric vector of damping factors to test (default: seq(0.5, 0.95, 0.05))
#' @param ltf_cutoff_range Numeric vector of LTF cutoffs to test (default: seq(0.80, 0.98, 0.02))
#' @param n_tfs_sample Number of TFs to sample for evaluation (default: 50, NULL for all)
#' @param top_n_targets Number of top targets per TF to use for evaluation (default: 100)
#' @param n_cores Number of cores for parallel processing (default: 1)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{best_params}{Named list with optimal damping_factor and ltf_cutoff}
#'   \item{best_score}{The consistency score for optimal parameters}
#'   \item{all_results}{Data frame with scores for all parameter combinations}
#'   \item{ligand_tf_coverage}{Data frame showing which TFs have upstream ligands}
#'
#' @details
#' The tuning process:
#' 1. Identifies TFs in the GRN that have upstream ligands (via lr_sig network)
#' 2. For each parameter combination, builds a ligand-target matrix
#' 3. Evaluates whether TF targets are enriched among ligand predictions
#' 4. Returns parameters that maximize this consistency
#'
#' This approach validates that signal correctly propagates:
#' Ligand -> (lr_sig) -> TF -> (GRN) -> Targets
#'
#' @examples
#' \dontrun{
#' grn <- import_grn("path/to/grn.csv")
#'
#' # Load NicheNet networks
#' lr_network <- readRDS("lr_network.rds")
#' weighted_networks <- readRDS("weighted_networks.rds")
#'
#' # Tune parameters
#' result <- tune_integration_params(
#'   grn = grn,
#'   lr_network = lr_network,
#'   weighted_networks = weighted_networks
#' )
#'
#' # Use optimal parameters
#' network <- integrate_networks(
#'   grn,
#'   damping_factor = result$best_params$damping_factor,
#'   ltf_cutoff = result$best_params$ltf_cutoff
#' )
#' }
#'
#' @export
tune_integration_params <- function(grn,
                                     lr_network,
                                     weighted_networks,
                                     damping_range = seq(0.5, 0.95, 0.05),
                                     ltf_cutoff_range = seq(0.80, 0.98, 0.02),
                                     n_tfs_sample = 50,
                                     top_n_targets = 100,
                                     n_cores = 1,
                                     verbose = TRUE) {

    # ── Validate inputs ──
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg nichenetr} is required for parameter tuning",
            "i" = "Install with: {.code devtools::install_github('saeyslab/nichenetr')}"
        ))
    }

    # Extract GRN edges
    if (inherits(grn, "ContextGRN")) {
        grn_edges <- grn@edges[, c("from", "to", "weight")]
    } else if (is.data.frame(grn)) {
        if (!all(c("from", "to", "weight") %in% colnames(grn))) {
            cli::cli_abort("GRN data.frame must have columns: from, to, weight")
        }
        grn_edges <- grn[, c("from", "to", "weight")]
    } else {
        cli::cli_abort("{.arg grn} must be a ContextGRN object or data.frame")
    }

    .validate_nichenet_networks(lr_network, weighted_networks)

    if (verbose) cli::cli_h2("Parameter Tuning for Network Integration")

    # ── Step 1: Find TFs with upstream ligands ──
    if (verbose) cli::cli_alert("Finding TFs with upstream ligands...")

    ligand_tf_map <- .find_ligand_tf_connections(
        lr_sig = weighted_networks$lr_sig,
        lr_network = lr_network,
        grn_tfs = unique(grn_edges$from),
        verbose = verbose
    )

    if (nrow(ligand_tf_map) == 0) {
        cli::cli_abort(c(
            "No TFs in the GRN have upstream ligands in the signaling network",
            "i" = "Check that TF names in GRN match gene names in NicheNet networks"
        ))
    }

    # Get TFs that have upstream ligands
    tfs_with_ligands <- unique(ligand_tf_map$tf)

    if (verbose) {
        cli::cli_alert_success("Found {.val {length(tfs_with_ligands)}} TFs with upstream ligands")
    }

    # Sample TFs if requested
    if (!is.null(n_tfs_sample) && length(tfs_with_ligands) > n_tfs_sample) {
        set.seed(42)  # For reproducibility
        tfs_with_ligands <- sample(tfs_with_ligands, n_tfs_sample)
        if (verbose) {
            cli::cli_alert_info("Sampled {.val {n_tfs_sample}} TFs for evaluation")
        }
    }

    # ── Step 2: Prepare TF-target ground truth ──
    if (verbose) cli::cli_alert("Preparing TF-target ground truth from GRN...")

    tf_targets_list <- .prepare_tf_targets(
        grn_edges = grn_edges,
        tfs = tfs_with_ligands,
        top_n = top_n_targets
    )

    # ── Step 3: Grid search over parameters ──
    param_grid <- expand.grid(
        damping_factor = damping_range,
        ltf_cutoff = ltf_cutoff_range,
        stringsAsFactors = FALSE
    )

    if (verbose) {
        cli::cli_alert("Testing {.val {nrow(param_grid)}} parameter combinations...")
        cli::cli_alert_info("Damping range: [{min(damping_range)}, {max(damping_range)}]")
        cli::cli_alert_info("LTF cutoff range: [{min(ltf_cutoff_range)}, {max(ltf_cutoff_range)}]")
    }

    # Prepare ligands list
    ligands <- unique(lr_network$from)
    ligands <- as.list(ligands)

    # Prepare weighted networks with custom GRN
    weighted_networks_custom <- list(
        lr_sig = weighted_networks$lr_sig,
        gr = grn_edges
    )

    # Progress bar
    if (verbose) {
        pb <- cli::cli_progress_bar(
            "Evaluating parameters",
            total = nrow(param_grid),
            clear = FALSE
        )
    }

    # Evaluate each parameter combination
    scores <- numeric(nrow(param_grid))

    for (i in seq_len(nrow(param_grid))) {
        d <- param_grid$damping_factor[i]
        l <- param_grid$ltf_cutoff[i]

        # Build ligand-target matrix
        ltm <- tryCatch({
            nichenetr::construct_ligand_target_matrix(
                weighted_networks = weighted_networks_custom,
                ligands = ligands,
                lr_network = lr_network,
                algorithm = "PPR",
                damping_factor = d,
                ltf_cutoff = l
            )
        }, error = function(e) {
            if (verbose) cli::cli_alert_warning("Error at d={d}, l={l}: {e$message}")
            return(NULL)
        })

        if (is.null(ltm)) {
            scores[i] <- NA
        } else {
            # Evaluate consistency
            scores[i] <- .evaluate_grn_consistency(
                ligand_target_matrix = as.matrix(ltm),
                ligand_tf_map = ligand_tf_map,
                tf_targets_list = tf_targets_list,
                tfs = tfs_with_ligands
            )
        }

        if (verbose) cli::cli_progress_update(id = pb)
    }

    if (verbose) cli::cli_progress_done(id = pb)

    # ── Step 4: Find best parameters ──
    param_grid$consistency_score <- scores

    best_idx <- which.max(param_grid$consistency_score)
    best_params <- list(
        damping_factor = param_grid$damping_factor[best_idx],
        ltf_cutoff = param_grid$ltf_cutoff[best_idx]
    )
    best_score <- param_grid$consistency_score[best_idx]

    if (verbose) {
        cli::cli_h3("Results")
        cli::cli_alert_success("Optimal parameters:")
        cli::cli_dl(c(
            "damping_factor" = "{.val {round(best_params$damping_factor, 3)}}",
            "ltf_cutoff" = "{.val {round(best_params$ltf_cutoff, 3)}}",
            "consistency_score" = "{.val {round(best_score, 4)}}"
        ))

        # Compare to NicheNet defaults
        nn_defaults <- c(damping_factor = 0.79, ltf_cutoff = 0.93)
        default_score <- param_grid$consistency_score[
            which.min(abs(param_grid$damping_factor - 0.79) +
                      abs(param_grid$ltf_cutoff - 0.93))
        ]

        cli::cli_alert_info("NicheNet defaults score: {.val {round(default_score, 4)}}")

        if (best_score > default_score) {
            improvement <- (best_score - default_score) / default_score * 100
            cli::cli_alert_success("Improvement over defaults: {.val {round(improvement, 1)}}%")
        }
    }

    list(
        best_params = best_params,
        best_score = best_score,
        all_results = param_grid,
        ligand_tf_coverage = ligand_tf_map
    )
}


#' Find ligand-TF connections through signaling network
#'
#' @param lr_sig Ligand-receptor-signaling network (from weighted_networks)
#' @param lr_network Ligand-receptor network
#' @param grn_tfs Character vector of TFs in the GRN
#' @param verbose Print messages
#'
#' @return Data frame with columns: ligand, tf
#'
#' @keywords internal
.find_ligand_tf_connections <- function(lr_sig, lr_network, grn_tfs, verbose = TRUE) {

    # Get ligands (sources in lr_network)
    ligands <- unique(lr_network$from)

    # The lr_sig network contains edges from ligands through signaling to TFs
    # We need to find which ligands can reach which TFs

    # Strategy: For each TF, check if it appears as a 'to' node that is

    # reachable from any ligand through the lr_sig network

    # Since lr_sig is large (~4M edges), we'll use an efficient approach:
    # 1. Filter lr_sig to edges where 'to' is a GRN TF
    # 2. Trace back to find ligands

    # Get all nodes that directly connect to TFs
    tf_incoming <- lr_sig[lr_sig$to %in% grn_tfs, ]

    if (nrow(tf_incoming) == 0) {
        return(data.frame(ligand = character(0), tf = character(0)))
    }

    # Check which of these 'from' nodes are ligands
    direct_ligand_tf <- tf_incoming[tf_incoming$from %in% ligands, c("from", "to")]
    colnames(direct_ligand_tf) <- c("ligand", "tf")

    # For indirect connections (ligand -> intermediate -> TF),
    # we trace one more step back
    intermediates <- setdiff(unique(tf_incoming$from), ligands)

    if (length(intermediates) > 0) {
        # Find ligands that connect to these intermediates
        intermediate_incoming <- lr_sig[lr_sig$to %in% intermediates &
                                         lr_sig$from %in% ligands, ]

        if (nrow(intermediate_incoming) > 0) {
            # Map: ligand -> intermediate -> tf
            # Join intermediate_incoming with tf_incoming
            indirect <- merge(
                intermediate_incoming[, c("from", "to")],
                tf_incoming[, c("from", "to")],
                by.x = "to",
                by.y = "from"
            )

            if (nrow(indirect) > 0) {
                indirect_ligand_tf <- data.frame(
                    ligand = indirect$from,
                    tf = indirect$to.y
                )
                direct_ligand_tf <- rbind(direct_ligand_tf, indirect_ligand_tf)
            }
        }
    }

    # Remove duplicates
    ligand_tf_map <- unique(direct_ligand_tf)

    if (verbose && nrow(ligand_tf_map) > 0) {
        n_ligands <- length(unique(ligand_tf_map$ligand))
        n_tfs <- length(unique(ligand_tf_map$tf))
        cli::cli_alert_info(
            "Found {.val {nrow(ligand_tf_map)}} ligand-TF connections ({.val {n_ligands}} ligands, {.val {n_tfs}} TFs)"
        )
    }

    ligand_tf_map
}


#' Prepare TF-target lists from GRN
#'
#' @param grn_edges GRN edges data.frame
#' @param tfs Character vector of TFs to include
#' @param top_n Number of top targets per TF (by weight)
#'
#' @return Named list: TF -> character vector of targets
#'
#' @keywords internal
.prepare_tf_targets <- function(grn_edges, tfs, top_n = 100) {

    tf_targets <- list()

    for (tf in tfs) {
        # Get targets for this TF, sorted by weight
        tf_edges <- grn_edges[grn_edges$from == tf, ]
        tf_edges <- tf_edges[order(-tf_edges$weight), ]

        # Take top N targets
        targets <- head(tf_edges$to, top_n)

        if (length(targets) > 0) {
            tf_targets[[tf]] <- targets
        }
    }

    tf_targets
}


#' Evaluate GRN propagation consistency
#'
#' @param ligand_target_matrix Ligand-target matrix
#' @param ligand_tf_map Data frame with ligand-TF connections
#' @param tf_targets_list Named list of TF -> targets
#' @param tfs TFs to evaluate
#'
#' @return Mean AUROC score
#'
#' @keywords internal
.evaluate_grn_consistency <- function(ligand_target_matrix,
                                       ligand_tf_map,
                                       tf_targets_list,
                                       tfs) {

    auroc_scores <- c()

    for (tf in tfs) {
        # Skip if TF has no targets
        if (!tf %in% names(tf_targets_list)) next

        targets <- tf_targets_list[[tf]]

        # Get ligands that signal through this TF
        upstream_ligands <- ligand_tf_map$ligand[ligand_tf_map$tf == tf]
        upstream_ligands <- intersect(upstream_ligands, colnames(ligand_target_matrix))

        if (length(upstream_ligands) == 0) next

        # For each upstream ligand, evaluate target prediction
        for (ligand in upstream_ligands) {
            # Get ligand's predicted scores for all genes
            if (!ligand %in% colnames(ligand_target_matrix)) next

            scores <- ligand_target_matrix[, ligand]

            # Filter to genes present in both
            common_genes <- intersect(names(scores), c(targets, names(scores)))
            if (length(common_genes) < 10) next

            scores <- scores[common_genes]

            # Create labels: TRUE if gene is a target of this TF
            labels <- common_genes %in% targets

            # Skip if no positive or negative examples
            if (sum(labels) == 0 || sum(!labels) == 0) next

            # Calculate AUROC
            auroc <- .calculate_auroc(scores, labels)

            if (!is.na(auroc)) {
                auroc_scores <- c(auroc_scores, auroc)
            }
        }
    }

    if (length(auroc_scores) == 0) {
        return(NA)
    }

    mean(auroc_scores, na.rm = TRUE)
}


#' Calculate AUROC
#'
#' @param scores Numeric vector of prediction scores
#' @param labels Logical vector of true labels
#'
#' @return AUROC value
#'
#' @keywords internal
.calculate_auroc <- function(scores, labels) {
    # Simple AUROC calculation using Mann-Whitney U statistic
    # AUROC = P(score of random positive > score of random negative)

    pos_scores <- scores[labels]
    neg_scores <- scores[!labels]

    n_pos <- length(pos_scores)
    n_neg <- length(neg_scores)

    if (n_pos == 0 || n_neg == 0) {
        return(NA)
    }

    # Count pairs where positive > negative
    # Use outer product for efficiency (but watch memory for large vectors)
    if (n_pos * n_neg > 1e7) {
        # Sample for large comparisons
        set.seed(123)
        sample_size <- min(1000, n_pos, n_neg)
        pos_sample <- sample(pos_scores, sample_size)
        neg_sample <- sample(neg_scores, sample_size)

        n_greater <- sum(outer(pos_sample, neg_sample, ">"))
        n_equal <- sum(outer(pos_sample, neg_sample, "=="))

        auroc <- (n_greater + 0.5 * n_equal) / (sample_size^2)
    } else {
        n_greater <- sum(outer(pos_scores, neg_scores, ">"))
        n_equal <- sum(outer(pos_scores, neg_scores, "=="))

        auroc <- (n_greater + 0.5 * n_equal) / (n_pos * n_neg)
    }

    auroc
}


#' Plot parameter tuning results
#'
#' @param tuning_result Result from tune_integration_params()
#' @param show_optimum Logical, highlight the optimal parameters (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' result <- tune_integration_params(grn, lr_network, weighted_networks)
#' plot_param_tuning(result)
#' }
#'
#' @export
plot_param_tuning <- function(tuning_result, show_optimum = TRUE) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        cli::cli_abort("ggplot2 is required for plotting")
    }

    df <- tuning_result$all_results
    best <- tuning_result$best_params

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = damping_factor, y = ltf_cutoff, fill = consistency_score)
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(
            name = "Consistency\nScore",
            option = "viridis",
            na.value = "grey50"
        ) +
        ggplot2::labs(
            title = "Parameter Tuning: GRN Propagation Consistency",
            x = "Damping Factor",
            y = "LTF Cutoff"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank()
        )

    if (show_optimum) {
        p <- p +
            ggplot2::annotate(
                "point",
                x = best$damping_factor,
                y = best$ltf_cutoff,
                shape = 4,
                size = 5,
                color = "red",
                stroke = 2
            ) +
            ggplot2::annotate(
                "text",
                x = best$damping_factor,
                y = best$ltf_cutoff,
                label = sprintf("Optimal\n(%.2f, %.2f)",
                               best$damping_factor, best$ltf_cutoff),
                vjust = -1,
                color = "red",
                size = 3
            )
    }

    p
}


#' Quick parameter recommendation based on GRN properties
#'
#' @description
#' Provides a quick parameter recommendation based on GRN network properties,
#' without running the full tuning procedure. Useful as a starting point.
#'
#' @param grn A ContextGRN object or data.frame with GRN edges
#' @param nichenet_gr NicheNet's gene regulatory network (for comparison)
#' @param verbose Logical, print recommendations
#'
#' @return Named list with recommended damping_factor and ltf_cutoff
#'
#' @examples
#' \dontrun{
#' grn <- import_grn("path/to/grn.csv")
#' weighted_networks <- readRDS("weighted_networks.rds")
#'
#' # Get quick recommendation
#' params <- recommend_integration_params(grn, weighted_networks$gr)
#' }
#'
#' @export
recommend_integration_params <- function(grn, nichenet_gr = NULL, verbose = TRUE) {

    # Extract GRN edges
    if (inherits(grn, "ContextGRN")) {
        grn_edges <- grn@edges
    } else if (is.data.frame(grn)) {
        grn_edges <- grn
    } else {
        cli::cli_abort("{.arg grn} must be a ContextGRN object or data.frame")
    }

    # NicheNet defaults
    base_damping <- 0.79
    base_ltf_cutoff <- 0.93

    if (is.null(nichenet_gr)) {
        if (verbose) {
            cli::cli_alert_info("No NicheNet GRN provided for comparison")
            cli::cli_alert("Returning NicheNet defaults")
        }
        return(list(
            damping_factor = base_damping,
            ltf_cutoff = base_ltf_cutoff,
            note = "Using NicheNet defaults (no comparison GRN provided)"
        ))
    }

    # Compare GRN properties
    if (verbose) cli::cli_h3("GRN Property Comparison")

    # Edge counts
    n_edges_custom <- nrow(grn_edges)
    n_edges_nn <- nrow(nichenet_gr)
    edge_ratio <- n_edges_custom / n_edges_nn

    # TF counts
    n_tfs_custom <- length(unique(grn_edges$from))
    n_tfs_nn <- length(unique(nichenet_gr$from))
    tf_ratio <- n_tfs_custom / n_tfs_nn

    # Target counts (unique targets)
    n_targets_custom <- length(unique(grn_edges$to))
    n_targets_nn <- length(unique(nichenet_gr$to))

    # Average degree (edges per TF)
    avg_degree_custom <- n_edges_custom / n_tfs_custom
    avg_degree_nn <- n_edges_nn / n_tfs_nn
    degree_ratio <- avg_degree_custom / avg_degree_nn

    if (verbose) {
        cli::cli_dl(c(
            "Custom GRN edges" = "{.val {format(n_edges_custom, big.mark = ',')}}",
            "NicheNet GRN edges" = "{.val {format(n_edges_nn, big.mark = ',')}}",
            "Edge ratio" = "{.val {round(edge_ratio, 2)}}",
            "Custom TFs" = "{.val {n_tfs_custom}}",
            "NicheNet TFs" = "{.val {n_tfs_nn}}",
            "Avg edges/TF (custom)" = "{.val {round(avg_degree_custom, 1)}}",
            "Avg edges/TF (NicheNet)" = "{.val {round(avg_degree_nn, 1)}}"
        ))
    }

    # Adjust parameters based on density
    # Denser network -> lower damping (signals spread more easily)
    # Sparser network -> higher damping (need more propagation)

    if (degree_ratio > 1.5) {
        # Much denser than NicheNet
        damping_adj <- -0.1
        ltf_adj <- 0.02
        reason <- "GRN is denser than NicheNet; reducing damping to limit over-propagation"
    } else if (degree_ratio < 0.5) {
        # Much sparser than NicheNet
        damping_adj <- 0.05
        ltf_adj <- -0.03
        reason <- "GRN is sparser than NicheNet; increasing damping for better coverage"
    } else {
        # Similar density
        damping_adj <- 0
        ltf_adj <- 0
        reason <- "GRN density similar to NicheNet; using default parameters"
    }

    recommended_damping <- min(0.95, max(0.5, base_damping + damping_adj))
    recommended_ltf <- min(0.98, max(0.80, base_ltf_cutoff + ltf_adj))

    if (verbose) {
        cli::cli_h3("Recommendation")
        cli::cli_alert_info(reason)
        cli::cli_dl(c(
            "Recommended damping_factor" = "{.val {round(recommended_damping, 2)}}",
            "Recommended ltf_cutoff" = "{.val {round(recommended_ltf, 2)}}"
        ))
        cli::cli_alert("For more precise tuning, use {.fn tune_integration_params}")
    }

    list(
        damping_factor = recommended_damping,
        ltf_cutoff = recommended_ltf,
        edge_ratio = edge_ratio,
        degree_ratio = degree_ratio,
        note = reason
    )
}

#' Show method for ContextBase objects
#'
#' @param object A ContextBase object
#' @return NULL (invisibly), prints to console
#'
#' @importFrom methods show
#' @keywords internal
setMethod("show", "ContextBase", function(object) {
    cli::cli_h1(class(object)[1])
    cli::cli_alert_info("Created: {format(object@provenance$created, '%Y-%m-%d %H:%M:%S')}")
    cli::cli_alert_info("Package version: {object@provenance$package_version}")
    if (!is.null(object@provenance$source)) {
        cli::cli_alert_info("Source: {object@provenance$source}")
    }
    invisible(NULL)
})


#' Show method for ContextGRN objects
#'
#' @param object A ContextGRN object
#' @return NULL (invisibly), prints to console
#'
#' @exportMethod show
setMethod("show", "ContextGRN", function(object) {
    n_edges <- nrow(object@edges)
    n_tfs <- length(unique(object@edges$from))
    n_targets <- length(unique(object@edges$to))

    cli::cli_h1("ContextGRN")
    cli::cli_text("")

    cli::cli_dl(c(
        "Edges" = "{.val {format(n_edges, big.mark = ',')}}",
        "TFs" = "{.val {n_tfs}}",
        "Targets" = "{.val {n_targets}}"
    ))

    # Show peak-gene links if present
    if (nrow(object@peak_gene_links) > 0) {
        cli::cli_dl(c(
            "Peak-gene links" = "{.val {format(nrow(object@peak_gene_links), big.mark = ',')}}"
        ))
    }

    # Show TF binding if present
    if (nrow(object@tf_binding) > 0) {
        cli::cli_dl(c(
            "TF binding sites" = "{.val {format(nrow(object@tf_binding), big.mark = ',')}}"
        ))
    }

    cli::cli_text("")
    cli::cli_alert_info("Source: {object@provenance$source}")
    cli::cli_alert_info("Created: {format(object@provenance$created, '%Y-%m-%d %H:%M:%S')}")

    # Show top TFs by edge count
    if (n_edges > 0) {
        tf_counts <- sort(table(object@edges$from), decreasing = TRUE)
        top_tfs <- head(names(tf_counts), 5)
        cli::cli_text("")
        cli::cli_alert_info("Top TFs (by target count): {.val {top_tfs}}")
    }

    invisible(NULL)
})


#' Show method for ContextNetwork objects
#'
#' @param object A ContextNetwork object
#' @return NULL (invisibly), prints to console
#'
#' @exportMethod show
setMethod("show", "ContextNetwork", function(object) {
    n_ligands <- ncol(object@ligand_target_matrix)
    n_targets <- nrow(object@ligand_target_matrix)

    cli::cli_h1("ContextNetwork")
    cli::cli_text("")

    cli::cli_dl(c(
        "Ligands" = "{.val {format(n_ligands, big.mark = ',')}}",
        "Targets" = "{.val {format(n_targets, big.mark = ',')}}"
    ))

    # Show integration info
    mode <- object@integration_params$mode %||% "unknown"
    cli::cli_dl(c(
        "Integration mode" = "{.val {mode}}"
    ))

    # Show weights if present
    if (!is.null(object@integration_params$weights)) {
        weights <- object@integration_params$weights
        weight_str <- paste(names(weights), round(weights, 2), sep = "=", collapse = ", ")
        cli::cli_dl(c(
            "Weights" = "{weight_str}"
        ))
    }

    # Show signaling network size
    if (nrow(object@signaling_network) > 0) {
        cli::cli_dl(c(
            "Signaling edges" = "{.val {format(nrow(object@signaling_network), big.mark = ',')}}"
        ))
    }

    # Show source GRN info
    if (!is.null(object@source_grn)) {
        cli::cli_text("")
        cli::cli_alert_info("Source GRN: {nrow(object@source_grn@edges)} edges")
    }

    cli::cli_text("")
    cli::cli_alert_info("Created: {format(object@provenance$created, '%Y-%m-%d %H:%M:%S')}")

    invisible(NULL)
})


#' Show method for ContextResult objects
#'
#' @param object A ContextResult object
#' @return NULL (invisibly), prints to console
#'
#' @exportMethod show
setMethod("show", "ContextResult", function(object) {
    cli::cli_h1("ContextResult")
    cli::cli_text("")

    # Analysis type
    type_label <- if (object@analysis_type == "multi_condition") {
        "Multi-condition (MultiNicheNet)"
    } else {
        "Standard (NicheNet)"
    }
    cli::cli_dl(c(
        "Analysis type" = "{type_label}"
    ))

    # Sender/receiver info
    if (!is.null(object@parameters$sender)) {
        cli::cli_dl(c(
            "Sender" = "{.val {object@parameters$sender}}",
            "Receiver" = "{.val {object@parameters$receiver}}"
        ))
    }

    # Condition info for multi-condition
    if (object@analysis_type == "multi_condition") {
        cli::cli_dl(c(
            "Condition" = "{.val {object@parameters$condition_oi}} vs {.val {object@parameters$condition_ref}}"
        ))
    }

    cli::cli_text("")

    # Results summary
    cli::cli_dl(c(
        "Top ligands" = "{.val {length(object@top_ligands)}}",
        "Top targets" = "{.val {length(object@top_targets)}}",
        "Ligand-target links" = "{.val {nrow(object@ligand_target_links)}}"
    ))

    # Show top ligands
    if (length(object@top_ligands) > 0) {
        cli::cli_text("")
        top5 <- head(object@top_ligands, 5)
        cli::cli_alert_info("Top 5 ligands: {.val {top5}}")

        # Show top ligand activity
        if (nrow(object@ligand_activities) > 0) {
            top_activity <- object@ligand_activities$aupr_corrected[1]
            cli::cli_alert_info("Top ligand AUPR: {round(top_activity, 4)}")
        }
    }

    cli::cli_text("")
    cli::cli_alert_info("Created: {format(object@provenance$created, '%Y-%m-%d %H:%M:%S')}")
    cli::cli_alert_info("Use {.fn plot} to visualize results")

    invisible(NULL)
})


#' Coerce ContextGRN to data.frame
#'
#' @param x A ContextGRN object
#' @return A data.frame of edges
#'
#' @export
setMethod("as.data.frame", "ContextGRN", function(x) {
    x@edges
})


#' Coerce ContextResult to data.frame
#'
#' @param x A ContextResult object
#' @return A data.frame of ligand activities
#'
#' @export
setMethod("as.data.frame", "ContextResult", function(x) {
    x@ligand_activities
})


#' Subset ContextGRN by TFs
#'
#' @param x A ContextGRN object
#' @param i TF names to subset
#' @param j Not used
#' @param ... Not used
#' @param drop Not used
#' @return A ContextGRN object with subset edges
#'
#' @export
setMethod("[", c("ContextGRN", "character"), function(x, i, j, ..., drop = TRUE) {
    edges <- x@edges[x@edges$from %in% i, , drop = FALSE]

    if (nrow(edges) == 0) {
        cli::cli_warn("No edges found for TFs: {.val {i}}")
    }

    ContextGRN(
        edges = edges,
        peak_gene_links = x@peak_gene_links,
        tf_binding = x@tf_binding,
        metadata = x@metadata
    )
})

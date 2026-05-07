#' Network Integration Functions
#'
#' @description
#' Functions for integrating custom GRNs with NicheNet prior knowledge networks
#' to create context-specific ligand-target matrices.
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#'
#' @name network-integrate
NULL


#' Integrate custom GRN with NicheNet networks
#'
#' @description Combine a context-specific GRN with NicheNet's prior knowledge
#'   networks to create an integrated ligand-target matrix.
#'
#' @param grn A ContextGRN object from import_grn()
#' @param lr_network Ligand-receptor network. If NULL, loads from NicheNet.
#' @param weighted_networks NicheNet weighted networks. If NULL, loads from NicheNet.
#' @param nichenet_database Path to NicheNet database directory (optional).
#'   If NULL, attempts to load from nichenetr package.
#' @param mode Integration mode: "replace" (default), "augment", or "weighted"
#' @param weights Named vector of weights for "weighted" mode (default: c(grn=0.5, nichenet=0.5))
#' @param damping_factor Damping factor for PPR algorithm (default: from nichenetr hyperparameters)
#' @param ltf_cutoff LTF cutoff for matrix construction (default: from nichenetr hyperparameters)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A ContextNetwork object containing the integrated ligand-target matrix
#'
#' @details
#' The integration process:
#' 1. Loads NicheNet's lr_network (ligand-receptor) and weighted_networks
#' 2. Replaces the gene regulatory (gr) component with the custom GRN
#' 3. Uses NicheNet's construct_ligand_target_matrix() with PPR algorithm
#'
#' Integration modes:
#' - "replace": Completely replace NicheNet's gr network with custom GRN
#' - "augment": Add custom GRN edges to NicheNet's gr network
#' - "weighted": Weighted combination of custom GRN and NicheNet gr
#'
#' @examples
#' \dontrun{
#' grn <- import_grn("path/to/gene_x_tf.csv")
#' network <- integrate_networks(grn)
#'
#' # With custom NicheNet database path
#' network <- integrate_networks(grn, nichenet_database = "/path/to/nichenet/")
#' }
#'
#' @export
integrate_networks <- function(grn,
                                lr_network = NULL,
                                weighted_networks = NULL,
                                nichenet_database = NULL,
                                mode = c("replace", "augment", "weighted"),
                                weights = c(grn = 0.5, nichenet = 0.5),
                                damping_factor = NULL,
                                ltf_cutoff = NULL,
                                verbose = TRUE) {
    # ── Validate inputs ──
    mode <- match.arg(mode)

    if (!inherits(grn, "ContextGRN")) {
        cli::cli_abort(c(
            "{.arg grn} must be a ContextGRN object",
            "i" = "Create one with {.fn import_grn}"
        ))
    }

    if (nrow(grn@edges) == 0) {
        cli::cli_abort("GRN has no edges")
    }

    # ── Check nichenetr dependency ──
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg nichenetr} is required for network integration",
            "i" = "Install with: {.code devtools::install_github('saeyslab/nichenetr')}"
        ))
    }

    # ── Load NicheNet networks ──
    if (verbose) cli::cli_h2("Network Integration")

    if (is.null(lr_network) || is.null(weighted_networks)) {
        if (verbose) cli::cli_alert("Loading NicheNet networks...")
        nn_data <- .load_nichenet_networks(nichenet_database, verbose)
        lr_network <- lr_network %||% nn_data$lr_network
        weighted_networks <- weighted_networks %||% nn_data$weighted_networks
    }

    # ── Validate NicheNet networks ──
    .validate_nichenet_networks(lr_network, weighted_networks)

    # ── Get ligands from lr_network ──
    # Logic from notebook: ligands <- lr_network %>% pull(from) %>% unique()
    ligands <- unique(lr_network$from)
    ligands <- as.list(ligands)

    if (verbose) {
        cli::cli_alert_info("Ligands from lr_network: {.val {length(ligands)}}")
    }

    # ── Prepare TF-target edges from GRN ──
    tf_target <- grn@edges[, c("from", "to", "weight")]

    if (verbose) {
        cli::cli_alert_info("GRN edges: {.val {format(nrow(tf_target), big.mark = ',')}}")
        cli::cli_alert_info("GRN TFs: {.val {length(unique(tf_target$from))}}")
        cli::cli_alert_info("GRN targets: {.val {length(unique(tf_target$to))}}")
    }

    # ── Create integrated weighted networks ──
    # Logic from notebook: weighted_networks_c <- list(lr_sig = ..., gr = tf_target)
    if (verbose) cli::cli_alert("Integrating networks (mode: {.val {mode}})...")

    weighted_networks_integrated <- .integrate_weighted_networks(
        weighted_networks = weighted_networks,
        tf_target = tf_target,
        mode = mode,
        weights = weights,
        verbose = verbose
    )

    # ── Get hyperparameters ──
    # Logic from notebook: uses hyperparameter_list from nichenetr
    hyperparameter_list <- nichenetr::hyperparameter_list

    if (is.null(damping_factor)) {
        damping_factor <- hyperparameter_list %>%
            dplyr::filter(.data$parameter == "damping_factor") %>%
            dplyr::pull(.data$avg_weight)
    }

    if (is.null(ltf_cutoff)) {
        ltf_cutoff <- hyperparameter_list %>%
            dplyr::filter(.data$parameter == "ltf_cutoff") %>%
            dplyr::pull(.data$avg_weight)
    }

    if (verbose) {
        cli::cli_alert_info("PPR damping_factor: {.val {round(damping_factor, 4)}}")
        cli::cli_alert_info("LTF cutoff: {.val {round(ltf_cutoff, 4)}}")
    }

    # ── Construct ligand-target matrix ──
    # Logic from notebook: construct_ligand_target_matrix() with PPR algorithm
    if (verbose) cli::cli_alert("Constructing ligand-target matrix (this may take a while)...")

    # nichenetr >= 2.2.1 added lr_network parameter; older versions don't have it
    ctm_args <- list(
        weighted_networks = weighted_networks_integrated,
        ligands = ligands,
        algorithm = "PPR",
        damping_factor = damping_factor,
        ltf_cutoff = ltf_cutoff
    )
    if ("lr_network" %in% names(formals(nichenetr::construct_ligand_target_matrix))) {
        ctm_args$lr_network <- lr_network
    }
    ligand_target_matrix <- do.call(nichenetr::construct_ligand_target_matrix, ctm_args)

    # ── Post-process matrix ──
    # Logic from notebook: fix column names and convert to matrix
    colnames(ligand_target_matrix) <- gsub("\\.", "-", colnames(ligand_target_matrix))
    ligand_target_matrix <- as.matrix(ligand_target_matrix)

    # ── Report results ──
    n_ligands <- ncol(ligand_target_matrix)
    n_targets <- nrow(ligand_target_matrix)

    if (verbose) {
        cli::cli_alert_success("Constructed ligand-target matrix:")
        cli::cli_dl(c(
            "Ligands" = "{.val {n_ligands}}",
            "Targets" = "{.val {n_targets}}",
            "Non-zero entries" = "{.val {format(sum(ligand_target_matrix > 0), big.mark = ',')}}"
        ))
    }

    # ── Create ContextNetwork object ──
    network <- new("ContextNetwork",
        ligand_target_matrix = ligand_target_matrix,
        ligand_receptor_matrix = matrix(nrow = 0, ncol = 0),
        signaling_network = weighted_networks_integrated$lr_sig,
        regulatory_network = weighted_networks_integrated$gr,
        weighted_networks = weighted_networks_integrated,
        integration_params = list(
            mode = mode,
            weights = if (mode == "weighted") weights else NULL,
            damping_factor = damping_factor,
            ltf_cutoff = ltf_cutoff,
            algorithm = "PPR"
        ),
        source_grn = grn,
        metadata = list(),
        provenance = .init_provenance(source = "integrate_networks")
    )

    network
}


#' Load NicheNet prior knowledge data
#'
#' @description Download NicheNet prior knowledge networks directly from Zenodo,
#'   matching the approach used in the official NicheNet and MultiNicheNet
#'   tutorials. Downloads are cached to \code{cache_dir} for reuse.
#'
#' @param cache_dir Directory to cache downloaded files. If NULL, files are
#'   downloaded each time. If a directory path is given, files are saved there
#'   and loaded from cache on subsequent calls.
#' @param version Which NicheNet data version to use:
#'   \describe{
#'     \item{"v2"}{Standard NicheNet v2 networks (Zenodo record 7074291).
#'       Includes lr_network, weighted_networks, and ligand_target_matrix.}
#'     \item{"v2_allinfo"}{Extended lr_network with additional metadata
#'       (Zenodo record 10229222). Used by MultiNicheNet tutorials.
#'       Only downloads lr_network; weighted_networks and ligand_target_matrix
#'       fall back to v2.}
#'   }
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A named list with:
#'   \describe{
#'     \item{lr_network}{Ligand-receptor network data.frame}
#'     \item{weighted_networks}{List with lr_sig and gr components}
#'     \item{ligand_target_matrix}{Ligand-target regulatory potential matrix}
#'   }
#'
#' @details
#' This function provides the same data that the NicheNet tutorials download
#' manually. Use the returned objects directly with \code{\link{integrate_networks}},
#' \code{\link{run_nichenet}}, or \code{\link{run_multinichenet}}:
#'
#' \preformatted{
#' # Download NicheNet data (cached for reuse)
#' nichenet_data <- load_nichenet_data(cache_dir = "~/.nichenet_data")
#'
#' # Use with integrate_networks
#' network <- integrate_networks(grn,
#'     lr_network = nichenet_data$lr_network,
#'     weighted_networks = nichenet_data$weighted_networks)
#'
#' # Use with run_nichenet
#' result <- run_nichenet(seurat_obj, network,
#'     lr_network = nichenet_data$lr_network,
#'     sender = "Macrophage", receiver = "Fibroblast",
#'     geneset_oi = de_genes)
#' }
#'
#' @examples
#' \dontrun{
#' # Download and cache NicheNet v2 data
#' nichenet_data <- load_nichenet_data(cache_dir = "~/.nichenet_data")
#'
#' # Use extended lr_network (MultiNicheNet style)
#' nichenet_data <- load_nichenet_data(version = "v2_allinfo",
#'                                     cache_dir = "~/.nichenet_data")
#' }
#'
#' @export
load_nichenet_data <- function(cache_dir = NULL,
                               version = c("v2", "v2_allinfo"),
                               verbose = TRUE) {
    version <- match.arg(version)

    # Zenodo URLs
    zenodo_v2 <- "https://zenodo.org/record/7074291/files"
    zenodo_allinfo <- "https://zenodo.org/record/10229222/files"

    files <- list(
        v2 = list(
            lr_network = list(
                url = paste0(zenodo_v2, "/lr_network_human_21122021.rds"),
                file = "lr_network_human_21122021.rds"
            ),
            weighted_networks = list(
                url = paste0(zenodo_v2, "/weighted_networks_nsga2r_final.rds"),
                file = "weighted_networks_nsga2r_final.rds"
            ),
            ligand_target_matrix = list(
                url = paste0(zenodo_v2, "/ligand_target_matrix_nsga2r_final.rds"),
                file = "ligand_target_matrix_nsga2r_final.rds"
            )
        ),
        v2_allinfo = list(
            lr_network = list(
                url = paste0(zenodo_allinfo,
                             "/lr_network_human_allInfo_30112033.rds"),
                file = "lr_network_human_allInfo_30112033.rds"
            )
        )
    )

    # Create cache directory if specified
    if (!is.null(cache_dir)) {
        cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE)
            if (verbose) cli::cli_alert_info("Created cache directory: {.path {cache_dir}}")
        }
    }

    # Helper to load a single file (cache or download)
    load_file <- function(file_info, label) {
        # Check cache first
        if (!is.null(cache_dir)) {
            cache_path <- file.path(cache_dir, file_info$file)
            if (file.exists(cache_path)) {
                if (verbose) cli::cli_alert_success("Loaded {label} from cache: {.file {cache_path}}")
                return(readRDS(cache_path))
            }
        }

        # Download
        if (verbose) cli::cli_alert("Downloading {label} from Zenodo...")
        data <- tryCatch({
            readRDS(url(file_info$url))
        }, error = function(e) {
            cli::cli_abort(c(
                "Failed to download {label} from Zenodo",
                "x" = e$message,
                "i" = "URL: {.url {file_info$url}}",
                "i" = "Check your internet connection or download manually"
            ))
        })
        if (verbose) cli::cli_alert_success("Downloaded {label}")

        # Cache if directory specified
        if (!is.null(cache_dir)) {
            tryCatch({
                saveRDS(data, file.path(cache_dir, file_info$file))
                if (verbose) cli::cli_alert_info("Cached to {.file {file.path(cache_dir, file_info$file)}}")
            }, error = function(e) {
                cli::cli_alert_warning("Could not cache file: {e$message}")
            })
        }

        data
    }

    result <- list()

    # Load lr_network (version-specific)
    if (version == "v2_allinfo") {
        result$lr_network <- load_file(files$v2_allinfo$lr_network, "lr_network (allInfo)")
    } else {
        result$lr_network <- load_file(files$v2$lr_network, "lr_network")
    }

    # weighted_networks and ligand_target_matrix always from v2
    result$weighted_networks <- load_file(files$v2$weighted_networks, "weighted_networks")
    result$ligand_target_matrix <- load_file(files$v2$ligand_target_matrix, "ligand_target_matrix")

    if (verbose) {
        cli::cli_alert_success("Loaded all NicheNet data ({version})")
        cli::cli_alert_info("lr_network: {nrow(result$lr_network)} interactions")
        if (is.matrix(result$ligand_target_matrix)) {
            cli::cli_alert_info(
                "ligand_target_matrix: {nrow(result$ligand_target_matrix)} ligands x {ncol(result$ligand_target_matrix)} targets"
            )
        }
    }

    result
}


#' Load NicheNet networks from database, Zenodo, or package
#'
#' @param nichenet_database Path to NicheNet database directory (optional)
#' @param verbose Logical, print messages
#'
#' @return List with lr_network and weighted_networks
#'
#' @keywords internal
.load_nichenet_networks <- function(nichenet_database = NULL, verbose = TRUE) {
    lr_network <- NULL
    weighted_networks <- NULL

    # Try loading from specified database path
    if (!is.null(nichenet_database)) {
        if (verbose) cli::cli_alert_info("Loading from: {.path {nichenet_database}}")

        # Try common file naming patterns
        lr_files <- c(
            "lr_network_human_21122021.rds",
            "lr_network_human.rds",
            "lr_network.rds"
        )

        wn_files <- c(
            "weighted_networks_nsga2r_final.rds",
            "weighted_networks.rds"
        )

        for (f in lr_files) {
            path <- file.path(nichenet_database, f)
            if (file.exists(path)) {
                lr_network <- readRDS(path)
                if (verbose) cli::cli_alert_success("Loaded lr_network from {.file {f}}")
                break
            }
        }

        for (f in wn_files) {
            path <- file.path(nichenet_database, f)
            if (file.exists(path)) {
                weighted_networks <- readRDS(path)
                if (verbose) cli::cli_alert_success("Loaded weighted_networks from {.file {f}}")
                break
            }
        }
    }

    # Try nichenetr package data
    if (is.null(lr_network)) {
        if (verbose) cli::cli_alert("Loading lr_network from nichenetr package...")
        lr_network <- tryCatch({
            nichenetr::lr_network
        }, error = function(e) {
            NULL
        })
    }

    if (is.null(weighted_networks)) {
        if (verbose) cli::cli_alert("Loading weighted_networks from nichenetr package...")
        weighted_networks <- tryCatch({
            nichenetr::weighted_networks
        }, error = function(e) {
            NULL
        })
    }

    # Fall back to Zenodo download
    if (is.null(lr_network) || is.null(weighted_networks)) {
        if (verbose) cli::cli_alert("Downloading missing NicheNet data from Zenodo...")
        zenodo <- tryCatch({
            load_nichenet_data(cache_dir = nichenet_database, verbose = verbose)
        }, error = function(e) {
            NULL
        })

        if (!is.null(zenodo)) {
            lr_network <- lr_network %||% zenodo$lr_network
            weighted_networks <- weighted_networks %||% zenodo$weighted_networks
        }
    }

    if (is.null(lr_network)) {
        cli::cli_abort(c(
            "Cannot load lr_network",
            "i" = "Use {.fn load_nichenet_data} to download, or provide path via {.arg nichenet_database}"
        ))
    }

    if (is.null(weighted_networks)) {
        cli::cli_abort(c(
            "Cannot load weighted_networks",
            "i" = "Use {.fn load_nichenet_data} to download, or provide path via {.arg nichenet_database}"
        ))
    }

    list(
        lr_network = lr_network,
        weighted_networks = weighted_networks
    )
}


#' Validate NicheNet networks structure
#'
#' @param lr_network Ligand-receptor network
#' @param weighted_networks Weighted networks list
#'
#' @keywords internal
.validate_nichenet_networks <- function(lr_network, weighted_networks) {
    # Validate lr_network
    if (!is.data.frame(lr_network)) {
        cli::cli_abort("lr_network must be a data.frame")
    }

    if (!"from" %in% colnames(lr_network)) {
        cli::cli_abort(c(
            "lr_network missing 'from' column",
            "i" = "Expected columns: from, to"
        ))
    }

    # Validate weighted_networks
    if (!is.list(weighted_networks)) {
        cli::cli_abort("weighted_networks must be a list")
    }

    required_components <- c("lr_sig", "gr")
    missing <- setdiff(required_components, names(weighted_networks))
    if (length(missing) > 0) {
        cli::cli_abort(c(
            "weighted_networks missing required components",
            "x" = "Missing: {.val {missing}}",
            "i" = "Expected: {.val {required_components}}"
        ))
    }

    invisible(TRUE)
}


#' Integrate weighted networks with custom GRN
#'
#' @param weighted_networks NicheNet weighted networks
#' @param tf_target Custom TF-target edges
#' @param mode Integration mode
#' @param weights Weights for weighted mode
#' @param verbose Print messages
#'
#' @return Integrated weighted networks list
#'
#' @keywords internal
.integrate_weighted_networks <- function(weighted_networks,
                                          tf_target,
                                          mode,
                                          weights,
                                          verbose) {
    # Logic from notebook: weighted_networks_c <- list(lr_sig = ..., gr = tf_target)

    if (mode == "replace") {
        # Replace gr component entirely with custom GRN
        weighted_networks_integrated <- list(
            lr_sig = weighted_networks$lr_sig,
            gr = tf_target
        )
    } else if (mode == "augment") {
        # Add custom GRN edges to existing gr network
        existing_gr <- weighted_networks$gr

        # Ensure column names match
        if (!all(c("from", "to", "weight") %in% colnames(existing_gr))) {
            cli::cli_abort("existing gr network must have from, to, weight columns")
        }

        # Combine and remove duplicates (keep custom GRN version)
        combined <- rbind(tf_target, existing_gr)
        combined <- combined[!duplicated(combined[, c("from", "to")]), ]

        weighted_networks_integrated <- list(
            lr_sig = weighted_networks$lr_sig,
            gr = combined
        )

        if (verbose) {
            cli::cli_alert_info("Augmented gr: {nrow(existing_gr)} + {nrow(tf_target)} -> {nrow(combined)} edges")
        }
    } else if (mode == "weighted") {
        # Weighted combination
        weights <- validate_weights(weights)

        existing_gr <- weighted_networks$gr

        # Merge on from/to, combine weights
        merged <- merge(
            tf_target,
            existing_gr,
            by = c("from", "to"),
            all = TRUE,
            suffixes = c(".grn", ".nichenet")
        )

        # Replace NA with 0
        merged$weight.grn[is.na(merged$weight.grn)] <- 0
        merged$weight.nichenet[is.na(merged$weight.nichenet)] <- 0

        # Weighted average
        merged$weight <- weights["grn"] * merged$weight.grn +
                         weights["nichenet"] * merged$weight.nichenet

        combined <- merged[, c("from", "to", "weight")]

        weighted_networks_integrated <- list(
            lr_sig = weighted_networks$lr_sig,
            gr = combined
        )

        if (verbose) {
            cli::cli_alert_info("Weighted combination: {nrow(combined)} edges")
        }
    }

    weighted_networks_integrated
}


#' Convert ContextNetwork to NicheNet input format
#'
#' @description Convert a ContextNetwork to the format expected by nichenetr
#'   functions. This allows users to follow standard NicheNet tutorials.
#'
#' @param network A ContextNetwork object
#'
#' @return A list containing:
#'   - ligand_target_matrix: Matrix of ligand-target regulatory potential
#'   - lr_network: Ligand-receptor network data.frame
#'   - sig_network: Signaling network data.frame
#'   - weighted_networks: List with lr_sig and gr components
#'
#' @examples
#' \dontrun{
#' grn <- import_grn("path/to/grn/")
#' network <- integrate_networks(grn)
#'
#' # Get NicheNet-compatible input
#' nn_input <- as_nichenet_input(network)
#'
#' # Use with nichenetr
#' library(nichenetr)
#' activities <- predict_ligand_activities(
#'     geneset = de_genes,
#'     background_expressed_genes = bg_genes,
#'     ligand_target_matrix = nn_input$ligand_target_matrix,
#'     potential_ligands = ligands
#' )
#' }
#'
#' @export
as_nichenet_input <- function(network) {
    if (!inherits(network, "ContextNetwork")) {
        cli::cli_abort(c(
            "{.arg network} must be a ContextNetwork object",
            "i" = "Create one with {.fn integrate_networks}"
        ))
    }

    list(
        ligand_target_matrix = network@ligand_target_matrix,
        lr_network = network@signaling_network,
        sig_network = network@signaling_network,
        weighted_networks = network@weighted_networks
    )
}

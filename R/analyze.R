#' High-Level Analysis Functions
#'
#' @description
#' User-friendly wrapper functions for cell-cell communication analysis that
#' return ContextResult objects. These wrap the lower-level run_nichenet() and
#' run_multinichenet() functions with additional validation, automatic DE
#' computation, and result packaging.
#'
#' After `integrate_networks()`, users have two options:
#' 1. Use these wrapper functions (recommended for most users)
#' 2. Follow NicheNet/MultiNicheNet pipelines directly with the integrated network
#'
#' @name analyze
NULL


#' Predict ligand activity using NicheNet framework
#'
#' @description Wrapper function for NicheNet ligand activity analysis. Identifies
#'   ligands from sender cells that regulate genes of interest in receiver cells.
#'
#'   This function follows the standard NicheNet workflow where differentially
#'   expressed genes are used as the target gene set. The DE comparison can be
#'   configured for two common scenarios:
#'
#'   1. **Condition comparison** (default): Compare receiver cells between conditions
#'      (e.g., disease vs control within fibroblasts)
#'   2. **Cell type comparison**: Compare different cell types directly
#'      (e.g., basal vs suprabasal epithelial cells)
#'
#' @param object A Seurat object containing expression data
#' @param network A ContextNetwork object from integrate_networks()
#' @param sender Character, name of sender cell type(s) in metadata. These are
#'   the cell types that express ligands of interest.
#' @param receiver Character, name of receiver cell type in metadata. This is
#'   the cell type where gene expression changes are measured.
#' @param geneset_oi Character vector of genes of interest in receiver cells.
#'   If NULL (default), DE genes are computed automatically. See Details for
#'   how DE genes are computed.
#' @param ident_col Column name in metadata for DE grouping (e.g., "condition",
#'   "cell_type"). REQUIRED when geneset_oi is NULL.
#' @param ident_1 Value(s) in ident_col for DE group 1 (upregulated genes).
#'   Can be a character vector to pool multiple values. REQUIRED when
#'   geneset_oi is NULL.
#' @param ident_2 Value(s) in ident_col for DE group 2 (reference).
#'   Can be a character vector to pool multiple values. REQUIRED when
#'   geneset_oi is NULL.
#' @param de_subset_receiver Logical, whether to subset to receiver cells before
#'   DE computation (default: TRUE).
#'   - TRUE: DE computed within receiver cells only (for condition comparisons)
#'   - FALSE: DE computed across all cells (for cell type comparisons)
#' @param group_col Column name in metadata containing cell type labels
#'   (default: "cell_type")
#' @param assay Assay to use from Seurat object (default: "RNA")
#' @param min_cells Minimum cells per cell type to include (default: 10)
#' @param min_pct Minimum percentage of cells expressing a gene to be considered
#'   expressed (default: 0.2, following NicheNet recommendation)
#' @param min_pct_de Minimum percentage of cells expressing a gene for
#'   FindMarkers DE computation (default: 0.05). Only used when computing DE
#'   genes automatically. Separate from min_pct to allow independent tuning.
#' @param logfc_threshold Log fold change threshold for DE gene selection
#'   (default: 0.25). Applied as a post-filter after FindMarkers (not passed
#'   to FindMarkers directly, so Bonferroni correction reflects all tested genes).
#' @param pval_threshold Adjusted p-value threshold for DE gene selection
#'   (default: 0.05). Only used when computing DE genes automatically.
#' @param top_n_ligands Number of top ligands to return (default: 30, following
#'   NicheNet recommendation)
#' @param top_n_targets Number of top targets per ligand (default: 100)
#' @param ligand_target_cutoff Minimum regulatory potential score to include a
#'   ligand-target link (default: 0). Set higher (e.g., 0.25) to retain only
#'   stronger predicted links.
#' @param lr_network Optional ligand-receptor network. If NULL, loads from
#'   nichenet_database or nichenetr package.
#' @param nichenet_database Path to NicheNet database directory containing
#'   lr_network_human.rds. Same path used in integrate_networks().
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A ContextResult object containing:
#'   \itemize{
#'     \item ligand_activities: data.frame of ligand activity scores (ranked by AUPR)
#'     \item ligand_target_links: data.frame of ligand-target relationships
#'     \item top_ligands: character vector of top-ranked ligands
#'     \item top_targets: character vector of top target genes
#'     \item receptor_info: data.frame of receptor information
#'     \item sender_receiver_info: list with sender/receiver cell type info
#'     \item expression_data: list of expression data used
#'     \item parameters: list of all parameters used
#'     \item analysis_type: "standard"
#'   }
#'
#' @details
#' ## Standard NicheNet Workflow
#'
#' This function implements the standard NicheNet workflow:
#'
#' 1. **Define sender and receiver cell types**
#' 2. **Identify DE genes** (between ident_1 and ident_2)
#' 3. **Identify expressed genes** in sender (potential ligands) and receiver
#' 4. **Run ligand activity prediction** (rank ligands by AUPR)
#' 5. **Extract top targets** regulated by top ligands
#'
#' ## DE Computation Modes
#'
#' The `de_subset_receiver` parameter controls how DE genes are computed:
#'
#' **de_subset_receiver = TRUE (default)**: For condition comparisons within a cell type
#' ```
#' # Example: Disease vs Control in Fibroblasts
#' receiver = "Fibroblast"
#' ident_col = "condition"
#' ident_1 = "Disease"
#' ident_2 = "Control"
#' # DE computed on Fibroblast cells only, comparing Disease vs Control
#' ```
#'
#' **de_subset_receiver = FALSE**: For cell type comparisons
#' ```
#' # Example: Basal vs Suprabasal epithelial cells
#' receiver = "basal"
#' ident_col = "cell_type"
#' ident_1 = "basal"
#' ident_2 = c("parabasal", "suprabasal")
#' de_subset_receiver = FALSE
#' # DE computed on ALL cells, comparing cell types directly
#' ```
#'
#' ## Using Your Own DE Genes
#'
#' To skip automatic DE computation, provide your own genes via `geneset_oi`.
#' This is useful when:
#' - You have custom DE results from a different method
#' - You want to focus on specific pathway genes
#' - You want to use marker genes instead of DE genes
#'
#' @examples
#' \dontrun{
#' # Import and integrate GRN
#' grn <- import_grn("path/to/grn/")
#' network <- integrate_networks(grn)
#'
#' # Example 1: Condition comparison (default behavior)
#' # Compare Disease vs Control WITHIN Fibroblast cells
#' result <- predict_ligand_activity_nichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = "Macrophage",
#'     receiver = "Fibroblast",
#'     ident_col = "disease_status",
#'     ident_1 = "Disease",
#'     ident_2 = "Control"
#'     # de_subset_receiver = TRUE (default)
#' )
#'
#' # Example 2: Cell type comparison
#' # Find basal cell markers vs parabasal/suprabasal
#' result <- predict_ligand_activity_nichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = all_celltypes,
#'     receiver = "basal",
#'     ident_col = "cell_type",
#'     ident_1 = "basal",
#'     ident_2 = c("parabasal", "suprabasal"),
#'     de_subset_receiver = FALSE  # Compare across all cells
#' )
#'
#' # Example 3: Multiple conditions pooled
#' result <- predict_ligand_activity_nichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = c("Macrophage", "Dendritic_cell"),
#'     receiver = "Fibroblast",
#'     ident_col = "disease_status",
#'     ident_1 = c("Disease_A", "Disease_B"),
#'     ident_2 = "Control"
#' )
#'
#' # Example 4: Custom gene set (skip DE computation)
#' pathway_genes <- c("IL6", "CCL2", "CXCL8", "MMP9")
#' result <- predict_ligand_activity_nichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = "Macrophage",
#'     receiver = "Fibroblast",
#'     geneset_oi = pathway_genes
#' )
#'
#' # Access results
#' top_ligands(result)
#' ligand_activities(result)
#' }
#'
#' @seealso
#' \code{\link{predict_ligand_activity_multinichenet}} for multi-condition analysis
#' \code{\link{run_nichenet}} for low-level NicheNet wrapper
#'
#' @export
predict_ligand_activity_nichenet <- function(object,
                                              network,
                                              sender,
                                              receiver,
                                              geneset_oi = NULL,
                                              ident_col = NULL,
                                              ident_1 = NULL,
                                              ident_2 = NULL,
                                              de_subset_receiver = TRUE,
                                              group_col = "cell_type",
                                              assay = "RNA",
                                              min_cells = 10,
                                              min_pct = 0.2,
                                              min_pct_de = 0.05,
                                              logfc_threshold = 0.25,
                                              pval_threshold = 0.05,
                                              top_n_ligands = 30,
                                              top_n_targets = 100,
                                              ligand_target_cutoff = 0,
                                              lr_network = NULL,
                                              nichenet_database = NULL,
                                              verbose = TRUE) {
    # ── Check dependencies ──
    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg nichenetr} is required for predict_ligand_activity_nichenet()",
            "i" = "Install with: {.code devtools::install_github('saeyslab/nichenetr')}"
        ))
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg Seurat} is required for predict_ligand_activity_nichenet()",
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
        cli::cli_abort("{.arg object} must be a Seurat object")
    }

    # Coerce sender/receiver to character (factors cause c() to return integers)
    sender <- as.character(sender)
    receiver <- as.character(receiver)

    # Validate group_col exists
    if (!group_col %in% colnames(object@meta.data)) {
        cli::cli_abort(c(
            "Column {.val {group_col}} not found in Seurat metadata",
            "i" = "Available columns: {.val {colnames(object@meta.data)}}"
        ))
    }

    # Validate cell types exist
    available_celltypes <- as.character(unique(object@meta.data[[group_col]]))
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

    # ── Validate DE parameters ──
    # If geneset_oi is NULL, require ident info for DE computation
    if (is.null(geneset_oi)) {
        if (is.null(ident_col) || is.null(ident_1) || is.null(ident_2)) {
            cli::cli_abort(c(
                "DE parameters required for automatic gene set computation",
                "",
                "Please provide ALL of:",
                "*" = "{.arg ident_col}: metadata column for DE grouping",
                "*" = "{.arg ident_1}: value(s) for DE group 1 (upregulated genes)",
                "*" = "{.arg ident_2}: value(s) for DE group 2 (reference)",
                "",
                "Or provide your own genes via {.arg geneset_oi} to skip DE computation"
            ))
        }
    }

    # Validate ident parameters if provided
    if (!is.null(ident_col)) {
        if (!ident_col %in% colnames(object@meta.data)) {
            cli::cli_abort(c(
                "Column {.val {ident_col}} not found in Seurat metadata",
                "i" = "Available columns: {.val {colnames(object@meta.data)}}"
            ))
        }

        available_idents <- unique(object@meta.data[[ident_col]])

        # ident_1 can be a vector
        if (!is.null(ident_1)) {
            missing_1 <- setdiff(ident_1, available_idents)
            if (length(missing_1) > 0) {
                cli::cli_abort(c(
                    "Value(s) not found in {.val {ident_col}}",
                    "x" = "Missing: {.val {missing_1}}",
                    "i" = "Available values: {.val {available_idents}}"
                ))
            }
        }

        # ident_2 can be a vector
        if (!is.null(ident_2)) {
            missing_2 <- setdiff(ident_2, available_idents)
            if (length(missing_2) > 0) {
                cli::cli_abort(c(
                    "Value(s) not found in {.val {ident_col}}",
                    "x" = "Missing: {.val {missing_2}}",
                    "i" = "Available values: {.val {available_idents}}"
                ))
            }
        }

        # Check no overlap
        if (!is.null(ident_1) && !is.null(ident_2)) {
            overlap <- intersect(ident_1, ident_2)
            if (length(overlap) > 0) {
                cli::cli_abort(c(
                    "Value(s) cannot be in both ident_1 and ident_2",
                    "x" = "Overlapping: {.val {overlap}}"
                ))
            }
        }
    }

    # Check minimum cells
    cell_counts <- table(object@meta.data[[group_col]])
    for (ct in c(sender, receiver)) {
        if (cell_counts[ct] < min_cells) {
            cli::cli_warn(c(
                "Cell type {.val {ct}} has fewer than {min_cells} cells",
                "i" = "Found {cell_counts[ct]} cells"
            ))
        }
    }

    if (verbose) {
        cli::cli_h1("NicheNet Ligand Activity Prediction")
        cli::cli_alert_info("Sender: {.val {paste(sender, collapse = ', ')}}")
        cli::cli_alert_info("Receiver: {.val {receiver}}")
        if (!is.null(ident_col)) {
            oi_str <- paste(ident_1, collapse = ", ")
            ref_str <- paste(ident_2, collapse = ", ")
            cli::cli_alert_info("DE comparison: {.val {oi_str}} vs {.val {ref_str}}")
            if (de_subset_receiver) {
                cli::cli_alert_info("DE mode: within {.val {receiver}} cells only")
            } else {
                cli::cli_alert_info("DE mode: across all cells")
            }
        }
    }

    # ── Compute DE genes if not provided ──
    if (is.null(geneset_oi)) {
        if (verbose) {
            oi_str <- paste(ident_1, collapse = ", ")
            ref_str <- paste(ident_2, collapse = ", ")
            if (de_subset_receiver) {
                cli::cli_alert("Computing DE genes: {receiver} cells, ({oi_str}) vs ({ref_str})")
            } else {
                cli::cli_alert("Computing DE genes: all cells, ({oi_str}) vs ({ref_str})")
            }
        }

        geneset_oi <- .compute_de_genes(
            object = object,
            receiver = receiver,
            group_col = group_col,
            ident_col = ident_col,
            ident_1 = ident_1,
            ident_2 = ident_2,
            de_subset_receiver = de_subset_receiver,
            min_pct = min_pct_de,
            logfc_threshold = logfc_threshold,
            pval_threshold = pval_threshold,
            verbose = verbose
        )

        if (length(geneset_oi) == 0) {
            oi_str <- paste(ident_1, collapse = ", ")
            ref_str <- paste(ident_2, collapse = ", ")
            cli::cli_abort(c(
                "No differentially expressed genes found",
                "i" = "Comparison: ({oi_str}) vs ({ref_str})",
                "",
                "Suggestions:",
                "*" = "Try lowering {.arg logfc_threshold} (current: {logfc_threshold})",
                "*" = "Try increasing {.arg pval_threshold} (current: {pval_threshold})",
                "*" = "Check that cells exist in both groups",
                "*" = "Provide your own genes via {.arg geneset_oi}"
            ))
        }

        if (verbose) {
            cli::cli_alert_success("Found {.val {length(geneset_oi)}} DE genes")
        }
    } else {
        if (!is.character(geneset_oi)) {
            cli::cli_abort("{.arg geneset_oi} must be a character vector of gene names")
        }
        if (verbose) {
            cli::cli_alert_info("Using provided geneset_oi: {.val {length(geneset_oi)}} genes")
        }
    }

    # ── Run NicheNet analysis ──
    nichenet_result <- run_nichenet(
        object = object,
        network = network,
        sender = sender,
        receiver = receiver,
        geneset_oi = geneset_oi,
        group_col = group_col,
        assay = assay,
        slot = "data",
        min_pct = min_pct,
        top_n_ligands = top_n_ligands,
        ligand_target_cutoff = ligand_target_cutoff,
        lr_network = lr_network,
        nichenet_database = nichenet_database,
        verbose = verbose
    )

    # ── Extract top targets ──
    if (verbose) cli::cli_alert("Extracting top targets...")

    ligand_target_matrix <- network@ligand_target_matrix
    top_ligands <- nichenet_result$top_ligands

    top_targets <- .get_top_targets(
        ligand_target_matrix = ligand_target_matrix,
        ligands = top_ligands,
        geneset_oi = geneset_oi,
        top_n = top_n_targets
    )

    # ── Create ContextResult ──
    if (verbose) cli::cli_alert("Creating ContextResult object...")

    result <- ContextResult(
        ligand_activities = nichenet_result$ligand_activities,
        ligand_target_links = nichenet_result$ligand_target_links,
        top_ligands = top_ligands,
        top_targets = top_targets,
        expression_data = list(
            sender_expressed = nichenet_result$nichenet_input$sender_expressed,
            receiver_expressed = nichenet_result$nichenet_input$receiver_expressed,
            background_genes = nichenet_result$nichenet_input$background_genes,
            geneset_oi = geneset_oi
        ),
        parameters = list(
            sender = sender,
            receiver = receiver,
            geneset_oi = geneset_oi,
            ident_col = ident_col,
            ident_1 = ident_1,
            ident_2 = ident_2,
            de_subset_receiver = de_subset_receiver,
            group_col = group_col,
            assay = assay,
            min_cells = min_cells,
            min_pct = min_pct,
            min_pct_de = min_pct_de,
            logfc_threshold = logfc_threshold,
            pval_threshold = pval_threshold,
            top_n_ligands = top_n_ligands,
            top_n_targets = top_n_targets,
            ligand_target_cutoff = ligand_target_cutoff
        ),
        analysis_type = "standard",
        receptor_info = nichenet_result$receptor_info,
        sender_receiver_info = list(
            sender = sender,
            receiver = receiver,
            sender_cell_count = sum(object@meta.data[[group_col]] %in% sender),
            receiver_cell_count = sum(object@meta.data[[group_col]] == receiver)
        )
    )

    if (verbose) {
        cli::cli_alert_success("Analysis complete")
        cli::cli_dl(c(
            "Top ligand" = "{.val {top_ligands[1]}}",
            "Ligands analyzed" = "{.val {nrow(nichenet_result$ligand_activities)}}",
            "Targets identified" = "{.val {length(top_targets)}}"
        ))
    }

    result
}


#' Predict ligand activity using MultiNicheNet framework
#'
#' @description Wrapper function for MultiNicheNet differential cell-cell
#'   communication analysis. Identifies ligands with differential activity
#'   between conditions, properly accounting for biological replicates.
#'
#'   MultiNicheNet performs sample-level analysis using pseudobulk aggregation,
#'   making it suitable for experiments with multiple samples per condition.
#'
#' @param object A Seurat object containing expression data with condition and
#'   sample information
#' @param network A ContextNetwork object from integrate_networks()
#' @param sender Character, name of sender cell type(s) in metadata. Can be a
#'   vector for multiple sender types.
#' @param receiver Character, name of receiver cell type in metadata. Analysis
#'   focuses on gene expression changes in this cell type.
#' @param ident_col Column name in metadata for DE grouping (e.g., "condition",
#'   "disease_status"). REQUIRED.
#' @param ident_1 Value in ident_col for condition of interest (e.g., "Disease").
#'   Genes upregulated in this group will be identified. REQUIRED.
#' @param ident_2 Value in ident_col for reference condition (e.g., "Control").
#'   REQUIRED.
#' @param sample_col Column name in metadata containing sample/replicate IDs
#'   (e.g., "patient_id", "sample"). REQUIRED for proper pseudobulk analysis.
#' @param geneset_oi Character vector of genes to focus analysis on. If NULL
#'   (default), DE genes are computed by MultiNicheNet per cell type.
#' @param group_col Column name in metadata containing cell type labels
#'   (default: "cell_type")
#' @param assay Assay to use from Seurat object (default: "RNA")
#' @param min_cells Minimum cells per sample/celltype/condition (default: 10)
#' @param top_n_ligands Number of top differential ligands to return (default: 30)
#' @param logfc_threshold Log fold change threshold for DE (default: 0.25).
#'   Passed to multinichenetr's `logFC_threshold`. Use 0.25 to match standard
#'   NicheNet/MultiNicheNet tutorials.
#' @param fraction_cutoff Minimum fraction of cells expressing a gene for it
#'   to be included in analysis (default: 0.05). Passed to multinichenetr.
#' @param top_n_target Top N target genes per ligand for activity scoring
#'   (default: 250). Passed to multinichenetr.
#' @param empirical_pval Logical, use permutation-based empirical p-values for
#'   DE (default: FALSE). Set to TRUE for permutation p-values; note this
#'   produces different ligand rankings than standard edgeR-based workflows.
#' @param lr_network Optional ligand-receptor network
#' @param nichenet_database Path to directory containing NicheNet data files
#' @param n_cores Number of cores for parallel processing (default: 4)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A ContextResult object containing:
#'   \itemize{
#'     \item ligand_activities: data.frame of ligand activity scores per condition
#'     \item ligand_target_links: data.frame of ligand-target relationships
#'     \item top_ligands: character vector of top differential ligands
#'     \item differential_results: list of DE results per cell type
#'     \item parameters: list of all parameters used
#'     \item analysis_type: "multi_condition"
#'   }
#'
#' @details
#' ## MultiNicheNet Workflow
#'
#' MultiNicheNet extends NicheNet for multi-sample, multi-condition experiments:
#'
#' 1. **Cell type filtering** - Identifies cell types with sufficient representation
#' 2. **Gene filtering** - Determines adequately expressed genes per cell type
#' 3. **Pseudobulk expression** - Calculates sample-level expression (accounts for replicates)
#' 4. **Differential expression** - Identifies DE genes between conditions per cell type
#' 5. **Ligand activity prediction** - Predicts ligand effects on receiver cells
#' 6. **Prioritization** - Ranks ligand-receptor pairs by multiple criteria
#'
#' ## Required Experimental Design
#'
#' MultiNicheNet requires:
#' - Multiple samples per condition (biological replicates)
#' - Cell type annotations
#' - Condition labels
#' - Sample/replicate identifiers
#'
#' If you have only one sample per condition, use `predict_ligand_activity_nichenet()`
#' instead.
#'
#' ## Contrast Definition
#'
#' The contrast compares `ident_1` vs `ident_2`. Positive differential
#' scores indicate higher activity in the condition of interest.
#'
#' @examples
#' \dontrun{
#' # Import and integrate GRN
#' grn <- import_grn("path/to/grn/")
#' network <- integrate_networks(grn)
#'
#' # Check your data structure
#' table(seurat_obj$condition, seurat_obj$sample)
#' #           Sample1 Sample2 Sample3 Sample4
#' # Control      1500    1400       0       0
#' # Disease         0       0    1600    1500
#'
#' # Run differential communication analysis
#' result <- predict_ligand_activity_multinichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = "Macrophage",
#'     receiver = "Fibroblast",
#'     ident_col = "condition",
#'     ident_1 = "Disease",
#'     ident_2 = "Control",
#'     sample_col = "sample"
#' )
#'
#' # Multiple sender cell types
#' result <- predict_ligand_activity_multinichenet(
#'     object = seurat_obj,
#'     network = network,
#'     sender = c("Macrophage", "Monocyte", "Dendritic_cell"),
#'     receiver = "Fibroblast",
#'     ident_col = "condition",
#'     ident_1 = "Disease",
#'     ident_2 = "Control",
#'     sample_col = "sample"
#' )
#'
#' # Access results
#' top_ligands(result)
#' ligand_activities(result)
#' }
#'
#' @seealso
#' \code{\link{predict_ligand_activity_nichenet}} for single-condition analysis
#' \code{\link{run_multinichenet}} for low-level MultiNicheNet wrapper
#'
#' @export
predict_ligand_activity_multinichenet <- function(object,
                                                   network,
                                                   sender,
                                                   receiver,
                                                   ident_col,
                                                   ident_1,
                                                   ident_2,
                                                   sample_col,
                                                   geneset_oi = NULL,
                                                   group_col = "cell_type",
                                                   assay = "RNA",
                                                   min_cells = 10,
                                                   top_n_ligands = 30,
                                                   logfc_threshold = 0.25,
                                                   fraction_cutoff = 0.05,
                                                   top_n_target = 250,
                                                   empirical_pval = FALSE,
                                                   lr_network = NULL,
                                                   nichenet_database = NULL,
                                                   n_cores = 4,
                                                   verbose = TRUE) {
    # ── Check dependencies ──
    if (!requireNamespace("multinichenetr", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg multinichenetr} is required for predict_ligand_activity_multinichenet()",
            "i" = "Install with: {.code devtools::install_github('saeyslab/multinichenetr')}"
        ))
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(c(
            "{.pkg Seurat} is required for predict_ligand_activity_multinichenet()",
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

    # Validate required condition parameters
    if (is.null(ident_col) || !nzchar(ident_col)) {
        cli::cli_abort(c(
            "{.arg ident_col} is required",
            "i" = "Specify the metadata column containing condition labels (e.g., 'condition')"
        ))
    }

    if (is.null(ident_1) || !nzchar(ident_1)) {
        cli::cli_abort(c(
            "{.arg ident_1} is required",
            "i" = "Specify the condition of interest (e.g., 'Disease')"
        ))
    }

    if (is.null(ident_2) || !nzchar(ident_2)) {
        cli::cli_abort(c(
            "{.arg ident_2} is required",
            "i" = "Specify the reference condition (e.g., 'Control')"
        ))
    }

    if (is.null(sample_col) || !nzchar(sample_col)) {
        cli::cli_abort(c(
            "{.arg sample_col} is required for MultiNicheNet",
            "i" = "Specify the metadata column containing sample/replicate IDs",
            "",
            "MultiNicheNet performs sample-level analysis to account for biological replicates.",
            "If you have only one sample per condition, use {.fn predict_ligand_activity_nichenet} instead."
        ))
    }

    if (!inherits(object, "Seurat")) {
        cli::cli_abort("{.arg object} must be a Seurat object")
    }

    # Coerce sender/receiver to character (factors cause c() to return integers)
    sender <- as.character(sender)
    receiver <- as.character(receiver)

    # Validate columns exist
    required_cols <- c(group_col, ident_col, sample_col)

    missing_cols <- setdiff(required_cols, colnames(object@meta.data))
    if (length(missing_cols) > 0) {
        cli::cli_abort(c(
            "Required columns not found in Seurat metadata",
            "x" = "Missing: {.val {missing_cols}}",
            "i" = "Available: {.val {colnames(object@meta.data)}}"
        ))
    }

    # Validate conditions exist
    available_conditions <- unique(object@meta.data[[ident_col]])
    if (!ident_1 %in% available_conditions) {
        cli::cli_abort(c(
            "Condition of interest not found in metadata",
            "x" = "{.val {ident_1}} not in {.val {ident_col}}",
            "i" = "Available conditions: {.val {available_conditions}}"
        ))
    }

    if (!ident_2 %in% available_conditions) {
        cli::cli_abort(c(
            "Reference condition not found in metadata",
            "x" = "{.val {ident_2}} not in {.val {ident_col}}",
            "i" = "Available conditions: {.val {available_conditions}}"
        ))
    }

    if (ident_1 == ident_2) {
        cli::cli_abort(c(
            "{.arg ident_1} and {.arg ident_2} must be different",
            "x" = "Both are {.val {ident_1}}"
        ))
    }

    # Validate cell types exist
    available_celltypes <- as.character(unique(object@meta.data[[group_col]]))
    all_celltypes <- unique(c(sender, receiver))
    missing_celltypes <- setdiff(all_celltypes, available_celltypes)
    if (length(missing_celltypes) > 0) {
        cli::cli_abort(c(
            "Cell types not found in metadata",
            "x" = "Missing: {.val {missing_celltypes}}",
            "i" = "Available: {.val {available_celltypes}}"
        ))
    }

    # Check minimum cells per cell type
    cell_counts <- table(object@meta.data[[group_col]])
    for (ct in all_celltypes) {
        if (!is.na(cell_counts[ct]) && cell_counts[ct] < min_cells) {
            cli::cli_warn(c(
                "Cell type {.val {ct}} has fewer than {min_cells} cells",
                "i" = "Found {cell_counts[ct]} cells"
            ))
        }
    }

    # Check samples per condition
    sample_condition <- unique(object@meta.data[, c(sample_col, ident_col)])
    samples_oi <- sum(sample_condition[[ident_col]] == ident_1)
    samples_ref <- sum(sample_condition[[ident_col]] == ident_2)

    if (samples_oi < 2 || samples_ref < 2) {
        cli::cli_warn(c(
            "MultiNicheNet works best with multiple samples per condition",
            "i" = "Samples in {.val {ident_1}}: {samples_oi}",
            "i" = "Samples in {.val {ident_2}}: {samples_ref}",
            "",
            "Consider using {.fn predict_ligand_activity_nichenet} if you have limited replicates."
        ))
    }

    if (verbose) {
        cli::cli_h1("MultiNicheNet Differential Communication Analysis")
        cli::cli_alert_info("Sender: {.val {paste(sender, collapse = ', ')}}")
        cli::cli_alert_info("Receiver: {.val {receiver}}")
        cli::cli_alert_info("Contrast: {.val {ident_1}} vs {.val {ident_2}}")
        cli::cli_alert_info("Samples: {samples_oi} ({ident_1}) vs {samples_ref} ({ident_2})")

        if (!is.null(geneset_oi)) {
            cli::cli_alert_info("Using provided geneset_oi: {.val {length(geneset_oi)}} genes")
        } else {
            cli::cli_alert_info("DE genes will be computed per cell type")
        }
    }

    if (!is.null(geneset_oi) && !is.character(geneset_oi)) {
        cli::cli_abort("{.arg geneset_oi} must be a character vector of gene names")
    }

    # ── Run MultiNicheNet analysis ──
    sender_receiver_pairs <- unlist(
        lapply(sender, function(s) lapply(receiver, function(r) c(s, r))),
        recursive = FALSE
    )

    multinichenet_result <- run_multinichenet(
        object = object,
        network = network,
        sender_receiver_pairs = sender_receiver_pairs,
        condition_col = ident_col,  # Translate to run_multinichenet parameter
        contrast = c(condition_oi = ident_1, condition_ref = ident_2),
        sample_col = sample_col,
        group_col = group_col,
        assay = assay,
        min_cells = min_cells,
        logfc_threshold = logfc_threshold,
        fraction_cutoff = fraction_cutoff,
        top_n_target = top_n_target,
        empirical_pval = empirical_pval,
        lr_network = lr_network,
        nichenet_database = nichenet_database,
        n_cores = n_cores,
        verbose = verbose
    )

    # ── Extract results ──
    # The result extraction mirrors the NicheNet pathway (predict_ligand_activity_nichenet)
    # in principle: both compute AUPR via nichenetr::predict_ligand_activities(),
    # extract ligand-target links, receptor info, and top targets. The key difference
    # is that MultiNicheNet ranks by prioritization_score (composite: AUPR + expression
    # + DE) rather than raw AUPR, and operates on pseudobulk per condition.
    if (verbose) cli::cli_alert("Processing results...")

    # ── Top ligands: from prioritization_table (ranked by prioritization_score) ──
    # NicheNet ranks by raw AUPR; MultiNicheNet ranks by prioritization_score which
    # combines AUPR (activity_scaled) with ligand/receptor expression and DE.
    top_ligands <- character()
    if (!is.null(multinichenet_result$prioritization_table) &&
        nrow(multinichenet_result$prioritization_table) > 0) {
        top_ligands <- utils::head(
            unique(multinichenet_result$prioritization_table$ligand),
            top_n_ligands
        )
    } else if (!is.null(multinichenet_result$ligand_activities) &&
               nrow(multinichenet_result$ligand_activities) > 0) {
        top_ligands <- utils::head(
            unique(multinichenet_result$ligand_activities$ligand),
            top_n_ligands
        )
    }

    # ── Geneset of interest: from receiver DE ──
    # NicheNet computes DE before analysis (Seurat FindMarkers);
    # MultiNicheNet computes pseudobulk DE internally (edgeR). We extract the
    # receiver DE genes post-hoc from the celltype_de table for top_targets.
    # celltype_de is a data.frame with columns: gene, cluster_id, logFC, p_val,
    # p_adj, contrast — filter by cluster_id to get receiver-specific DE genes.
    if (is.null(geneset_oi)) {
        geneset_oi <- character()
        celltype_de <- multinichenet_result$celltype_de
        if (!is.null(celltype_de) && is.data.frame(celltype_de) &&
            nrow(celltype_de) > 0) {
            # Filter to receiver cell types and condition-of-interest contrast
            contrast_oi <- paste0(ident_1, "-", ident_2)
            recv_de <- celltype_de[
                celltype_de$cluster_id %in% receiver &
                celltype_de$contrast == contrast_oi, ]
            if (nrow(recv_de) > 0) {
                geneset_oi <- unique(recv_de$gene[recv_de$p_adj < 0.05 &
                                                   recv_de$logFC > 0])
            }
        }
        if (length(geneset_oi) == 0 && verbose) {
            cli::cli_alert_warning(
                "No significant DE genes found in receiver cell types for top_targets"
            )
        }
    }

    # ── Top targets: same as NicheNet pathway ──
    ligand_target_matrix <- network@ligand_target_matrix
    top_targets <- character()
    if (length(top_ligands) > 0 && length(geneset_oi) > 0) {
        top_targets <- .get_top_targets(
            ligand_target_matrix = ligand_target_matrix,
            ligands = top_ligands,
            geneset_oi = geneset_oi,
            top_n = 100
        )
    }

    # ── Create ContextResult ──
    if (verbose) cli::cli_alert("Creating ContextResult object...")

    # ── Ligand activities: summarize per ligand ──
    # multinichenetr's ligand_activities has one row per (ligand, target, receiver,
    # contrast, direction_regulation). We add test_ligand and aupr_corrected aliases
    # for compatibility with NicheNet's format.
    ligand_activities <- multinichenet_result$ligand_activities
    if (is.null(ligand_activities)) {
        ligand_activities <- data.frame(
            test_ligand = character(),
            aupr_corrected = numeric()
        )
    } else {
        if (!"test_ligand" %in% colnames(ligand_activities)) {
            if ("ligand" %in% colnames(ligand_activities)) {
                ligand_activities$test_ligand <- ligand_activities$ligand
            }
        }
        if (!"aupr_corrected" %in% colnames(ligand_activities)) {
            if ("activity" %in% colnames(ligand_activities)) {
                ligand_activities$aupr_corrected <- ligand_activities$activity
            } else {
                ligand_activities$aupr_corrected <- 0
            }
        }
    }

    # ── Ligand-target links: extract from ligand_activities ──
    # In multinichenetr, target info is embedded in the ligand_activities table
    # (columns: target, ligand_target_weight) — unlike NicheNet which uses
    # get_active_ligand_target_df() separately. Extract and reshape to match
    # NicheNet format: ligand, target, weight.
    ligand_target_links <- data.frame()
    if (!is.null(ligand_activities) && nrow(ligand_activities) > 0 &&
        all(c("ligand", "target", "ligand_target_weight") %in%
            colnames(ligand_activities))) {
        # Filter to top ligands and condition of interest
        contrast_oi <- paste0(ident_1, "-", ident_2)
        ltl <- ligand_activities[
            ligand_activities$ligand %in% top_ligands &
            ligand_activities$contrast == contrast_oi &
            ligand_activities$direction_regulation == "up" &
            !is.na(ligand_activities$target) &
            !is.na(ligand_activities$ligand_target_weight), ]
        if (nrow(ltl) > 0) {
            ligand_target_links <- unique(ltl[, c("ligand", "target",
                                                   "ligand_target_weight")])
            colnames(ligand_target_links)[
                colnames(ligand_target_links) == "ligand_target_weight"] <- "weight"
        }
    }

    # Extract receptor_info from prioritization table
    receptor_info <- data.frame()
    if (!is.null(multinichenet_result$prioritization_table) &&
        nrow(multinichenet_result$prioritization_table) > 0) {
        prio <- multinichenet_result$prioritization_table
        if (all(c("ligand", "receptor") %in% colnames(prio))) {
            receptor_info <- unique(prio[, c("ligand", "receptor"), drop = FALSE])
            colnames(receptor_info) <- c("from", "to")
        }
    }

    result <- ContextResult(
        ligand_activities = ligand_activities,
        ligand_target_links = ligand_target_links,
        top_ligands = top_ligands,
        top_targets = top_targets,
        expression_data = list(
            geneset_oi = geneset_oi,
            group_prioritization_tbl =
                multinichenet_result$multinichenet_output$prioritization_tables$group_prioritization_tbl
        ),
        parameters = list(
            sender = sender,
            receiver = receiver,
            geneset_oi = geneset_oi,
            condition_oi = ident_1,
            ident_col = ident_col,
            ident_1 = ident_1,
            ident_2 = ident_2,
            sample_col = sample_col,
            group_col = group_col,
            assay = assay,
            min_cells = min_cells,
            top_n_ligands = top_n_ligands,
            logfc_threshold = logfc_threshold,
            fraction_cutoff = fraction_cutoff,
            top_n_target = top_n_target,
            empirical_pval = empirical_pval
        ),
        analysis_type = "multi_condition",
        receptor_info = receptor_info,
        sender_receiver_info = list(
            sender = sender,
            receiver = receiver,
            sender_receiver_pairs = sender_receiver_pairs,
            sender_cell_count = sum(object@meta.data[[group_col]] %in% sender),
            receiver_cell_count = sum(object@meta.data[[group_col]] %in% receiver)
        ),
        differential_results = multinichenet_result$celltype_de
    )

    if (verbose) {
        cli::cli_alert_success("Analysis complete")
        cli::cli_dl(c(
            "Top differential ligands" = "{.val {length(top_ligands)}}",
            "Contrast" = "{.val {ident_1}} vs {.val {ident_2}}"
        ))
    }

    result
}


#' Compute DE genes for NicheNet analysis
#'
#' @description Computes differentially expressed genes using Seurat's FindMarkers.
#'   Supports two modes:
#'   - de_subset_receiver = TRUE: DE computed within receiver cells only
#'   - de_subset_receiver = FALSE: DE computed across all cells
#'
#' @param object Seurat object
#' @param receiver Receiver cell type
#' @param group_col Cell type column
#' @param ident_col Column for DE grouping
#' @param ident_1 Value(s) for DE group 1 (upregulated genes)
#' @param ident_2 Value(s) for DE group 2 (reference)
#' @param de_subset_receiver Whether to subset to receiver cells first
#' @param min_pct Minimum fraction of cells expressing a gene for FindMarkers
#' @param logfc_threshold Log fold change threshold
#' @param pval_threshold P-value threshold
#' @param verbose Print progress
#'
#' @return Character vector of DE gene names
#'
#' @keywords internal
.compute_de_genes <- function(object,
                               receiver,
                               group_col,
                               ident_col,
                               ident_1,
                               ident_2,
                               de_subset_receiver = TRUE,
                               min_pct = 0.05,
                               logfc_threshold = 0.25,
                               pval_threshold = 0.05,
                               verbose = TRUE) {

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort("Seurat is required for DE computation")
    }

    # Subset based on de_subset_receiver mode
    if (de_subset_receiver) {
        # Mode 1: DE within receiver cells only (condition comparison)
        receiver_cells <- which(object@meta.data[[group_col]] == receiver)

        if (length(receiver_cells) == 0) {
            cli::cli_abort("No cells found for receiver cell type {.val {receiver}}")
        }

        object_sub <- object[, receiver_cells]
    } else {
        # Mode 2: DE across all cells (cell type comparison)
        object_sub <- object
    }

    # Check cells per ident (handle vectors)
    ident_counts <- table(object_sub@meta.data[[ident_col]])

    # Count cells in ident_1 (can be vector)
    cells_1 <- sum(ident_counts[names(ident_counts) %in% ident_1], na.rm = TRUE)
    if (cells_1 < 3) {
        location_str <- if (de_subset_receiver) paste0("in ", receiver) else "in dataset"
        cli::cli_abort(c(
            "Insufficient cells in ident_1 {location_str}",
            "x" = "Values: {.val {ident_1}}",
            "i" = "Found: {cells_1} cells total"
        ))
    }

    # Count cells in ident_2 (can be vector)
    cells_2 <- sum(ident_counts[names(ident_counts) %in% ident_2], na.rm = TRUE)
    if (cells_2 < 3) {
        location_str <- if (de_subset_receiver) paste0("in ", receiver) else "in dataset"
        cli::cli_abort(c(
            "Insufficient cells in ident_2 {location_str}",
            "x" = "Values: {.val {ident_2}}",
            "i" = "Found: {cells_2} cells total"
        ))
    }

    # Create a grouped column for FindMarkers when multiple values are provided.
    # Pool cells into "group_1" vs "group_2" groups.
    if (length(ident_1) > 1 || length(ident_2) > 1) {
        # Create grouping variable
        object_sub$de_group <- ifelse(
            object_sub@meta.data[[ident_col]] %in% ident_1,
            "group_1",
            ifelse(
                object_sub@meta.data[[ident_col]] %in% ident_2,
                "group_2",
                NA
            )
        )
        # Subset to cells in either group
        cells_in_comparison <- which(!is.na(object_sub$de_group))
        object_sub <- object_sub[, cells_in_comparison]
        Seurat::Idents(object_sub) <- object_sub$de_group
        de_ident_1 <- "group_1"
        de_ident_2 <- "group_2"
    } else {
        # Single value on each side
        Seurat::Idents(object_sub) <- object_sub@meta.data[[ident_col]]
        de_ident_1 <- ident_1
        de_ident_2 <- ident_2
    }

    # Increase future.globals.maxSize for large Seurat objects
    old_maxsize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = Inf)
    on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)

    # Run DE: ident_1 vs ident_2
    # Do NOT pass logfc.threshold to FindMarkers — use Seurat's default (0.1).
    # Passing a higher threshold (e.g. 0.25) reduces the number of genes tested,
    # which weakens Bonferroni correction and inflates geneset_oi. Post-filtering
    # by logfc_threshold is applied below instead, matching the standard NicheNet
    # notebook workflow.
    de_results <- Seurat::FindMarkers(
        object_sub,
        ident.1 = de_ident_1,
        ident.2 = de_ident_2,
        min.pct = min_pct,
        verbose = FALSE
    )

    # Filter DE genes
    if (nrow(de_results) == 0) {
        return(character())
    }

    # Use avg_log2FC if available, otherwise avg_logFC
    fc_col <- if ("avg_log2FC" %in% colnames(de_results)) "avg_log2FC" else "avg_logFC"
    pval_col <- if ("p_val_adj" %in% colnames(de_results)) "p_val_adj" else "p_val"

    # Filter by thresholds (upregulated genes in ident_1)
    sig_genes <- rownames(de_results)[
        de_results[[fc_col]] >= logfc_threshold &
        de_results[[pval_col]] < pval_threshold
    ]

    sig_genes
}


#' Get top targets for ligands
#'
#' @param ligand_target_matrix Ligand-target regulatory potential matrix
#' @param ligands Character vector of ligands
#' @param geneset_oi Genes of interest to filter targets
#' @param top_n Number of top targets to return
#'
#' @return Character vector of top target genes
#'
#' @keywords internal
.get_top_targets <- function(ligand_target_matrix,
                              ligands,
                              geneset_oi,
                              top_n = 100) {
    # Filter to ligands in matrix
    ligands_in_matrix <- intersect(ligands, colnames(ligand_target_matrix))

    if (length(ligands_in_matrix) == 0) {
        return(character())
    }

    # Filter to geneset_oi targets
    targets_in_matrix <- intersect(geneset_oi, rownames(ligand_target_matrix))

    if (length(targets_in_matrix) == 0) {
        return(character())
    }

    # Get submatrix
    submat <- ligand_target_matrix[targets_in_matrix, ligands_in_matrix, drop = FALSE]

    # Get max regulatory potential per target across ligands
    target_scores <- apply(submat, 1, max)

    # Sort and get top targets
    sorted_targets <- names(sort(target_scores, decreasing = TRUE))
    utils::head(sorted_targets, top_n)
}

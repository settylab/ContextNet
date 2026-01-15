#' @title ContextGRN Class
#' @description Container for gene regulatory network data. Stores edges (TF-target
#'   relationships), optional peak-gene links, and TF binding information.
#'
#' @slot edges A data.frame with columns: from (TF), to (target), weight
#' @slot nodes A data.frame with gene information
#' @slot regulons A list of regulons (TF -> target gene sets)
#' @slot peak_gene_links A data.frame linking peaks to genes
#' @slot tf_binding A data.frame with TF binding information
#' @slot metadata A list of user-defined metadata (inherited from ContextBase)
#' @slot provenance A list of provenance information (inherited from ContextBase)
#'
#' @name ContextGRN-class
#' @rdname ContextGRN-class
#' @exportClass ContextGRN
setClass("ContextGRN",
    contains = "ContextBase",
    slots = c(
        edges = "data.frame",
        nodes = "data.frame",
        regulons = "list",
        peak_gene_links = "data.frame",
        tf_binding = "data.frame"
    ),
    prototype = list(
        edges = data.frame(
            from = character(),
            to = character(),
            weight = numeric(),
            stringsAsFactors = FALSE
        ),
        nodes = data.frame(
            gene = character(),
            stringsAsFactors = FALSE
        ),
        regulons = list(),
        peak_gene_links = data.frame(),
        tf_binding = data.frame(),
        metadata = list(),
        provenance = list()
    )
)


#' Validity check for ContextGRN objects
#' @keywords internal
setValidity("ContextGRN", function(object) {
    errors <- character()

    # Check required columns in edges
    required_cols <- c("from", "to", "weight")
    if (!all(required_cols %in% colnames(object@edges))) {
        missing <- setdiff(required_cols, colnames(object@edges))
        errors <- c(errors, sprintf(
            "edges data.frame missing required columns: %s",
            paste(missing, collapse = ", ")
        ))
    }

    # Check weight range if edges exist
    if (nrow(object@edges) > 0 && "weight" %in% colnames(object@edges)) {
        weights <- object@edges$weight
        if (any(is.na(weights))) {
            errors <- c(errors, "edges$weight contains NA values")
        } else if (any(weights < 0 | weights > 1)) {
            errors <- c(errors, "edges$weight must be in range [0, 1]")
        }
    }

    # Check nodes if edges exist
    if (nrow(object@edges) > 0 && nrow(object@nodes) > 0) {
        all_genes <- unique(c(object@edges$from, object@edges$to))
        if (!all(all_genes %in% object@nodes$gene)) {
            errors <- c(errors, "Not all genes in edges are present in nodes")
        }
    }

    if (length(errors) == 0) TRUE else errors
})


#' Create a ContextGRN object
#'
#' @param edges A data.frame with columns: from (TF), to (target), weight.
#'   Column names are flexible and will be normalized.
#' @param peak_gene_links Optional data.frame linking peaks to genes
#' @param tf_binding Optional data.frame with TF binding information
#' @param metadata Optional list of user-defined metadata
#'
#' @return A ContextGRN object
#'
#' @examples
#' \dontrun{
#' edges <- data.frame(
#'     from = c("TF1", "TF1", "TF2"),
#'     to = c("Gene1", "Gene2", "Gene1"),
#'     weight = c(0.8, 0.6, 0.7)
#' )
#' grn <- ContextGRN(edges)
#' }
#'
#' @export
ContextGRN <- function(edges,
                       peak_gene_links = NULL,
                       tf_binding = NULL,
                       metadata = list()) {
    # Normalize edge column names
    edges <- .normalize_edge_colnames(edges)

    # Validate edges
    .validate_edges(edges)

    # Build nodes from edges
    nodes <- data.frame(
        gene = unique(c(edges$from, edges$to)),
        stringsAsFactors = FALSE
    )

    # Build regulons from edges
    regulons <- .build_regulons(edges)

    # Handle optional inputs
    if (is.null(peak_gene_links)) {
        peak_gene_links <- data.frame()
    }
    if (is.null(tf_binding)) {
        tf_binding <- data.frame()
    }

    # Create object
    new("ContextGRN",
        edges = edges,
        nodes = nodes,
        regulons = regulons,
        peak_gene_links = peak_gene_links,
        tf_binding = tf_binding,
        metadata = metadata,
        provenance = .init_provenance(source = "manual")
    )
}


#' Normalize edge data.frame column names
#'
#' @param edges A data.frame with edge information
#' @return A data.frame with standardized column names
#'
#' @keywords internal
.normalize_edge_colnames <- function(edges) {
    if (!is.data.frame(edges)) {
        cli::cli_abort("{.arg edges} must be a data.frame")
    }

    cols <- colnames(edges)
    cols_lower <- tolower(cols)

    # Map common column name variants
    from_patterns <- c("from", "tf", "source", "regulator", "transcription_factor")
    to_patterns <- c("to", "target", "gene", "dest", "destination")
    weight_patterns <- c("weight", "score", "confidence", "correlation", "cor")

    new_cols <- cols

    # Find 'from' column
    from_match <- which(cols_lower %in% from_patterns)
    if (length(from_match) == 0) {
        cli::cli_abort(c(
            "Cannot find 'from' column in edges",
            "i" = "Expected one of: {.val {from_patterns}}"
        ))
    }
    new_cols[from_match[1]] <- "from"

    # Find 'to' column
    to_match <- which(cols_lower %in% to_patterns)
    if (length(to_match) == 0) {
        cli::cli_abort(c(
            "Cannot find 'to' column in edges",
            "i" = "Expected one of: {.val {to_patterns}}"
        ))
    }
    new_cols[to_match[1]] <- "to"

    # Find 'weight' column
    weight_match <- which(cols_lower %in% weight_patterns)
    if (length(weight_match) == 0) {
        cli::cli_abort(c(
            "Cannot find 'weight' column in edges",
            "i" = "Expected one of: {.val {weight_patterns}}"
        ))
    }
    new_cols[weight_match[1]] <- "weight"

    colnames(edges) <- new_cols
    edges[, c("from", "to", "weight")]
}


#' Validate edge data.frame
#'
#' @param edges A data.frame with edge information
#' @return TRUE if valid, otherwise throws an error
#'
#' @keywords internal
.validate_edges <- function(edges) {
    if (nrow(edges) == 0) {
        cli::cli_abort("edges data.frame is empty")
    }

    if (any(is.na(edges$from) | is.na(edges$to))) {
        cli::cli_abort("edges contains NA values in 'from' or 'to' columns")
    }

    if (any(is.na(edges$weight))) {
        cli::cli_abort("edges contains NA values in 'weight' column")
    }

    if (any(edges$weight < 0 | edges$weight > 1)) {
        cli::cli_abort(c(
            "Edge weights must be in range [0, 1]",
            "x" = "Found weights in range [{min(edges$weight)}, {max(edges$weight)}]"
        ))
    }

    invisible(TRUE)
}


#' Build regulons from edges
#'
#' @param edges A data.frame with from, to, weight columns
#' @return A named list of regulons (TF -> character vector of targets)
#'
#' @keywords internal
.build_regulons <- function(edges) {
    tfs <- unique(edges$from)
    regulons <- lapply(tfs, function(tf) {
        edges$to[edges$from == tf]
    })
    names(regulons) <- tfs
    regulons
}


#' Get the number of edges in a ContextGRN
#'
#' @param x A ContextGRN object
#' @return The number of edges
#'
#' @export
setMethod("length", "ContextGRN", function(x) {
    nrow(x@edges)
})


#' Get TFs from a ContextGRN
#'
#' @param object A ContextGRN object
#' @return Character vector of TF names
#'
#' @export
#' @docType methods
#' @rdname tfs-methods
setGeneric("tfs", function(object) standardGeneric("tfs"))


#' @rdname tfs-methods
#' @aliases tfs,ContextGRN-method
setMethod("tfs", "ContextGRN", function(object) {
    unique(object@edges$from)
})


#' Get target genes from a ContextGRN
#'
#' @param object A ContextGRN object
#' @return Character vector of target gene names
#'
#' @export
#' @docType methods
#' @rdname targets-methods
setGeneric("targets", function(object) standardGeneric("targets"))


#' @rdname targets-methods
#' @aliases targets,ContextGRN-method
setMethod("targets", "ContextGRN", function(object) {
    unique(object@edges$to)
})


#' Get edges from a ContextGRN
#'
#' @param object A ContextGRN object
#' @return A data.frame with edge information
#'
#' @export
#' @docType methods
#' @rdname edges-methods
setGeneric("edges", function(object) standardGeneric("edges"))


#' @rdname edges-methods
#' @aliases edges,ContextGRN-method
setMethod("edges", "ContextGRN", function(object) {
    object@edges
})


#' Get regulons from a ContextGRN
#'
#' @param object A ContextGRN object
#' @return A named list of regulons
#'
#' @export
#' @docType methods
#' @rdname regulons-methods
setGeneric("regulons", function(object) standardGeneric("regulons"))


#' @rdname regulons-methods
#' @aliases regulons,ContextGRN-method
setMethod("regulons", "ContextGRN", function(object) {
    object@regulons
})

#' @title ContextNetwork Class
#' @description Container for integrated network combining custom GRN with
#'   NicheNet prior knowledge networks.
#'
#' @slot ligand_target_matrix A matrix of ligand-target regulatory potential
#' @slot ligand_receptor_matrix A matrix of ligand-receptor interactions
#' @slot signaling_network A data.frame of signaling interactions
#' @slot regulatory_network A data.frame of regulatory interactions
#' @slot weighted_networks A list of weighted network components
#' @slot integration_params A list of parameters used for integration
#' @slot source_grn The source ContextGRN object (or NULL)
#' @slot metadata A list of user-defined metadata (inherited from ContextBase)
#' @slot provenance A list of provenance information (inherited from ContextBase)
#'
#' @name ContextNetwork-class
#' @rdname ContextNetwork-class
#' @exportClass ContextNetwork
setClass("ContextNetwork",
    contains = "ContextBase",
    slots = c(
        ligand_target_matrix = "matrix",
        ligand_receptor_matrix = "matrix",
        signaling_network = "data.frame",
        regulatory_network = "data.frame",
        weighted_networks = "list",
        integration_params = "list",
        source_grn = "ANY"  # ContextGRN or NULL
    ),
    prototype = list(
        ligand_target_matrix = matrix(nrow = 0, ncol = 0),
        ligand_receptor_matrix = matrix(nrow = 0, ncol = 0),
        signaling_network = data.frame(),
        regulatory_network = data.frame(),
        weighted_networks = list(),
        integration_params = list(),
        source_grn = NULL,
        metadata = list(),
        provenance = list()
    )
)


#' Validity check for ContextNetwork objects
#' @keywords internal
setValidity("ContextNetwork", function(object) {
    errors <- character()

    # Check ligand_target_matrix
    if (!is.matrix(object@ligand_target_matrix)) {
        errors <- c(errors, "ligand_target_matrix must be a matrix")
    }

    # Check source_grn type if not NULL
    if (!is.null(object@source_grn) && !inherits(object@source_grn, "ContextGRN")) {
        errors <- c(errors, "source_grn must be a ContextGRN object or NULL")
    }

    # Check integration_params has required fields
    if (length(object@integration_params) > 0) {
        required_params <- c("mode")
        missing_params <- setdiff(required_params, names(object@integration_params))
        if (length(missing_params) > 0) {
            errors <- c(errors, sprintf(
                "integration_params missing required fields: %s",
                paste(missing_params, collapse = ", ")
            ))
        }
    }

    if (length(errors) == 0) TRUE else errors
})


#' Create a ContextNetwork object
#'
#' @param ligand_target_matrix A matrix of ligand-target regulatory potential
#' @param ligand_receptor_matrix A matrix of ligand-receptor interactions
#' @param signaling_network A data.frame of signaling interactions
#' @param regulatory_network A data.frame of regulatory interactions
#' @param weighted_networks A list of weighted network components
#' @param integration_params A list of parameters used for integration
#' @param source_grn Optional ContextGRN object that was used as source
#' @param metadata Optional list of user-defined metadata
#'
#' @return A ContextNetwork object
#'
#' @export
ContextNetwork <- function(ligand_target_matrix,
                           ligand_receptor_matrix = NULL,
                           signaling_network = NULL,
                           regulatory_network = NULL,
                           weighted_networks = NULL,
                           integration_params = list(mode = "weighted"),
                           source_grn = NULL,
                           metadata = list()) {
    # Validate ligand_target_matrix
    if (!is.matrix(ligand_target_matrix)) {
        cli::cli_abort("{.arg ligand_target_matrix} must be a matrix")
    }

    # Handle optional inputs
    if (is.null(ligand_receptor_matrix)) {
        ligand_receptor_matrix <- matrix(nrow = 0, ncol = 0)
    }
    if (is.null(signaling_network)) {
        signaling_network <- data.frame()
    }
    if (is.null(regulatory_network)) {
        regulatory_network <- data.frame()
    }
    if (is.null(weighted_networks)) {
        weighted_networks <- list()
    }

    # Ensure integration_params has mode
    if (!"mode" %in% names(integration_params)) {
        integration_params$mode <- "weighted"
    }

    # Create object
    new("ContextNetwork",
        ligand_target_matrix = ligand_target_matrix,
        ligand_receptor_matrix = ligand_receptor_matrix,
        signaling_network = signaling_network,
        regulatory_network = regulatory_network,
        weighted_networks = weighted_networks,
        integration_params = integration_params,
        source_grn = source_grn,
        metadata = metadata,
        provenance = .init_provenance(source = "integrate_networks")
    )
}


#' Get ligands from a ContextNetwork
#'
#' @param object A ContextNetwork object
#' @return Character vector of ligand names
#'
#' @export
#' @docType methods
#' @rdname ligands-methods
setGeneric("ligands", function(object) standardGeneric("ligands"))


#' @rdname ligands-methods
#' @aliases ligands,ContextNetwork-method
setMethod("ligands", "ContextNetwork", function(object) {
    colnames(object@ligand_target_matrix)
})


#' Get target genes from a ContextNetwork
#'
#' @rdname targets-methods
#' @aliases targets,ContextNetwork-method
setMethod("targets", "ContextNetwork", function(object) {
    rownames(object@ligand_target_matrix)
})


#' Get integration parameters from a ContextNetwork
#'
#' @param object A ContextNetwork object
#' @return A list of integration parameters
#'
#' @export
#' @docType methods
#' @rdname integration_params-methods
setGeneric("integration_params", function(object) standardGeneric("integration_params"))


#' @rdname integration_params-methods
#' @aliases integration_params,ContextNetwork-method
setMethod("integration_params", "ContextNetwork", function(object) {
    object@integration_params
})


#' Get the ligand-target matrix from a ContextNetwork
#'
#' @param object A ContextNetwork object
#' @return A matrix of ligand-target regulatory potential
#'
#' @export
#' @docType methods
#' @rdname ligand_target_matrix-methods
setGeneric("ligand_target_matrix", function(object) standardGeneric("ligand_target_matrix"))


#' @rdname ligand_target_matrix-methods
#' @aliases ligand_target_matrix,ContextNetwork-method
setMethod("ligand_target_matrix", "ContextNetwork", function(object) {
    object@ligand_target_matrix
})

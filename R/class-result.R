#' @title ContextResult Class
#' @description Container for cell-cell communication analysis results.
#'   Stores ligand activities, ligand-target links, and analysis parameters.
#'
#' @slot ligand_activities A data.frame of ligand activity scores
#' @slot ligand_target_links A data.frame of ligand-target relationships
#' @slot receptor_info A data.frame of receptor information
#' @slot sender_receiver_info A list with sender/receiver cell type info
#' @slot top_ligands Character vector of top-ranked ligands
#' @slot top_targets Character vector of top target genes
#' @slot expression_data A list containing expression data used in analysis
#' @slot differential_results Optional list of differential analysis results
#' @slot parameters A list of parameters used for analysis
#' @slot analysis_type Character string: "standard" or "multi_condition"
#' @slot metadata A list of user-defined metadata (inherited from ContextBase)
#' @slot provenance A list of provenance information (inherited from ContextBase)
#'
#' @name ContextResult-class
#' @rdname ContextResult-class
#' @exportClass ContextResult
setClass("ContextResult",
    contains = "ContextBase",
    slots = c(
        ligand_activities = "data.frame",
        ligand_target_links = "data.frame",
        receptor_info = "data.frame",
        sender_receiver_info = "list",
        top_ligands = "character",
        top_targets = "character",
        expression_data = "list",
        differential_results = "ANY",  # list or NULL
        parameters = "list",
        analysis_type = "character"
    ),
    prototype = list(
        ligand_activities = data.frame(),
        ligand_target_links = data.frame(),
        receptor_info = data.frame(),
        sender_receiver_info = list(),
        top_ligands = character(),
        top_targets = character(),
        expression_data = list(),
        differential_results = NULL,
        parameters = list(),
        analysis_type = "standard",
        metadata = list(),
        provenance = list()
    )
)


#' Validity check for ContextResult objects
#' @keywords internal
setValidity("ContextResult", function(object) {
    errors <- character()

    # Check analysis_type
    valid_types <- c("standard", "multi_condition")
    if (!object@analysis_type %in% valid_types) {
        errors <- c(errors, sprintf(
            "analysis_type must be one of: %s",
            paste(valid_types, collapse = ", ")
        ))
    }

    # Check ligand_activities has required columns
    if (nrow(object@ligand_activities) > 0) {
        required_cols <- c("test_ligand", "aupr_corrected")
        missing_cols <- setdiff(required_cols, colnames(object@ligand_activities))
        if (length(missing_cols) > 0) {
            errors <- c(errors, sprintf(
                "ligand_activities missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            ))
        }
    }

    if (length(errors) == 0) TRUE else errors
})


#' Create a ContextResult object
#'
#' @param ligand_activities A data.frame of ligand activity scores
#' @param ligand_target_links A data.frame of ligand-target relationships
#' @param top_ligands Character vector of top-ranked ligands
#' @param top_targets Character vector of top target genes
#' @param expression_data A list containing expression data used in analysis
#' @param parameters A list of parameters used for analysis
#' @param analysis_type Character: "standard" or "multi_condition"
#' @param receptor_info Optional data.frame of receptor information
#' @param sender_receiver_info Optional list with sender/receiver info
#' @param differential_results Optional list of differential analysis results
#' @param metadata Optional list of user-defined metadata
#'
#' @return A ContextResult object
#'
#' @export
ContextResult <- function(ligand_activities,
                          ligand_target_links,
                          top_ligands,
                          top_targets,
                          expression_data,
                          parameters,
                          analysis_type = "standard",
                          receptor_info = NULL,
                          sender_receiver_info = NULL,
                          differential_results = NULL,
                          metadata = list()) {
    # Handle optional inputs
    if (is.null(receptor_info)) {
        receptor_info <- data.frame()
    }
    if (is.null(sender_receiver_info)) {
        sender_receiver_info <- list()
    }

    # Determine source based on analysis type
    source <- if (analysis_type == "multi_condition") {
        "analyze_communication_multi"
    } else {
        "analyze_communication"
    }

    # Create object
    new("ContextResult",
        ligand_activities = ligand_activities,
        ligand_target_links = ligand_target_links,
        receptor_info = receptor_info,
        sender_receiver_info = sender_receiver_info,
        top_ligands = top_ligands,
        top_targets = top_targets,
        expression_data = expression_data,
        differential_results = differential_results,
        parameters = parameters,
        analysis_type = analysis_type,
        metadata = metadata,
        provenance = .init_provenance(source = source)
    )
}


#' Get ligand activities from a ContextResult
#'
#' @param object A ContextResult object
#' @return A data.frame of ligand activity scores
#'
#' @export
#' @docType methods
#' @rdname ligand_activities-methods
setGeneric("ligand_activities", function(object) standardGeneric("ligand_activities"))


#' @rdname ligand_activities-methods
#' @aliases ligand_activities,ContextResult-method
setMethod("ligand_activities", "ContextResult", function(object) {
    object@ligand_activities
})


#' Get top ligands from a ContextResult
#'
#' @param object A ContextResult object
#' @param n Number of top ligands to return (default: all)
#' @return Character vector of top ligand names
#'
#' @export
#' @docType methods
#' @rdname top_ligands-methods
setGeneric("top_ligands", function(object, n = NULL) standardGeneric("top_ligands"))


#' @rdname top_ligands-methods
#' @aliases top_ligands,ContextResult-method
setMethod("top_ligands", "ContextResult", function(object, n = NULL) {
    ligands <- object@top_ligands
    if (!is.null(n) && n < length(ligands)) {
        ligands <- ligands[seq_len(n)]
    }
    ligands
})


#' Get top targets from a ContextResult
#'
#' @param object A ContextResult object
#' @param n Number of top targets to return (default: all)
#' @return Character vector of top target gene names
#'
#' @export
#' @docType methods
#' @rdname top_targets-methods
setGeneric("top_targets", function(object, n = NULL) standardGeneric("top_targets"))


#' @rdname top_targets-methods
#' @aliases top_targets,ContextResult-method
setMethod("top_targets", "ContextResult", function(object, n = NULL) {
    targets <- object@top_targets
    if (!is.null(n) && n < length(targets)) {
        targets <- targets[seq_len(n)]
    }
    targets
})


#' Get ligand-target links from a ContextResult
#'
#' @param object A ContextResult object
#' @return A data.frame of ligand-target links
#'
#' @export
#' @docType methods
#' @rdname ligand_target_links-methods
setGeneric("ligand_target_links", function(object) standardGeneric("ligand_target_links"))


#' @rdname ligand_target_links-methods
#' @aliases ligand_target_links,ContextResult-method
setMethod("ligand_target_links", "ContextResult", function(object) {
    object@ligand_target_links
})


#' Get analysis parameters from a ContextResult
#'
#' @param object A ContextResult object
#' @return A list of analysis parameters
#'
#' @export
#' @docType methods
#' @rdname parameters-methods
setGeneric("parameters", function(object) standardGeneric("parameters"))


#' @rdname parameters-methods
#' @aliases parameters,ContextResult-method
setMethod("parameters", "ContextResult", function(object) {
    object@parameters
})

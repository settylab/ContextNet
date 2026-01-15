#' @title ContextBase Class
#' @description Virtual base class for all ContextNet objects. Provides common
#'   slots and methods for metadata and provenance tracking.
#'
#' @slot metadata A list containing user-defined metadata
#' @slot provenance A list containing creation time, package version, and source info
#'
#' @name ContextBase-class
#' @rdname ContextBase-class
#' @exportClass ContextBase
setClass("ContextBase",
    contains = "VIRTUAL",
    slots = c(
        metadata = "list",
        provenance = "list"
    ),
    prototype = list(
        metadata = list(),
        provenance = list(
            created = NULL,
            package_version = NULL
        )
    )
)


#' Get provenance information from a ContextNet object
#'
#' @param object A ContextBase object
#' @return A list containing provenance information
#'
#' @export
#' @docType methods
#' @rdname provenance-methods
setGeneric("provenance", function(object) standardGeneric("provenance"))


#' @rdname provenance-methods
#' @aliases provenance,ContextBase-method
setMethod("provenance", "ContextBase", function(object) {
    object@provenance
})


#' Get metadata from a ContextNet object
#'
#' @param object A ContextBase object
#' @return A list containing metadata
#'
#' @export
#' @docType methods
#' @rdname metadata-methods
setGeneric("metadata", function(object) standardGeneric("metadata"))


#' @rdname metadata-methods
#' @aliases metadata,ContextBase-method
setMethod("metadata", "ContextBase", function(object) {
    object@metadata
})


#' Set metadata for a ContextNet object
#'
#' @param object A ContextBase object
#' @param value A list of metadata to set
#' @return The modified object
#'
#' @export
#' @docType methods
#' @rdname metadata-methods
setGeneric("metadata<-", function(object, value) standardGeneric("metadata<-"))


#' @rdname metadata-methods
#' @aliases metadata<-,ContextBase-method
setMethod("metadata<-", "ContextBase", function(object, value) {
    if (!is.list(value)) {
        cli::cli_abort("{.arg value} must be a list")
    }
    object@metadata <- value
    object
})


#' Initialize provenance for a new object
#'
#' @param source Character string describing the source of the object
#' @param source_version Optional version string for the source
#' @return A list with provenance information
#'
#' @keywords internal
.init_provenance <- function(source = "manual", source_version = NA_character_) {
    pkg_version <- tryCatch(
        as.character(utils::packageVersion("ContextNet")),
        error = function(e) "0.0.0.9000"
    )

    list(
        created = Sys.time(),
        package_version = pkg_version,
        source = source,
        source_version = source_version
    )
}

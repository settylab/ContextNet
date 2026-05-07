#' AnnData to Seurat Conversion
#'
#' @description
#' Convert AnnData (.h5ad) files to Seurat objects using the native R
#' anndataR package (no Python required).
#'
#' @name interop-anndata
NULL


#' Convert AnnData file to Seurat object
#'
#' @description
#' Read an h5ad file and convert it to a Seurat object. Uses the anndataR
#' package for native R reading without Python dependencies.
#'
#' @param path Path to the h5ad file
#' @param assay_name Name for the Seurat assay (default: "RNA")
#'
#' @return A Seurat object
#'
#' @examples
#' \dontrun{
#' # Convert AnnData to Seurat
#' seurat_obj <- anndata_to_seurat("data/adata.h5ad")
#'
#' # With custom assay name
#' seurat_obj <- anndata_to_seurat("data/adata.h5ad", assay_name = "SCT")
#' }
#'
#' @export
anndata_to_seurat <- function(path, assay_name = "RNA") {
    # Validate path
    validate_path(path, must_exist = TRUE, extension = "h5ad")

    # Check for required packages
    if (!requireNamespace("anndataR", quietly = TRUE)) {
        cli::cli_abort(c(
            "Package {.pkg anndataR} is required for h5ad conversion",
            "i" = "Install with: {.code install.packages('anndataR')}"
        ))
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(c(
            "Package {.pkg Seurat} is required for conversion",
            "i" = "Install with: {.code install.packages('Seurat')}"
        ))
    }

    cli::cli_alert_info("Reading h5ad file: {.path {path}}")

    # Read AnnData using anndataR
    adata <- tryCatch({
        anndataR::read_h5ad(path)
    }, error = function(e) {
        cli::cli_abort(c(
            "Failed to read h5ad file",
            "x" = conditionMessage(e),
            "i" = "Check that the file is a valid h5ad format"
        ))
    })

    cli::cli_alert_info("Converting to Seurat object")

    # Convert to Seurat
    seurat_obj <- tryCatch({
        anndataR::to_Seurat(adata, assay_name = assay_name)
    }, error = function(e) {
        cli::cli_abort(c(
            "Failed to convert AnnData to Seurat",
            "x" = conditionMessage(e)
        ))
    })

    cli::cli_alert_info("Conversion complete: {ncol(seurat_obj)} cells, {nrow(seurat_obj)} features")

    seurat_obj
}

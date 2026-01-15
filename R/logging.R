#' Logging and Verbosity System
#'
#' @description
#' Configurable logging system for ContextNet. Supports different log levels
#' and optional file output.
#'
#' @name logging
NULL


#' Setup ContextNet logging
#'
#' @description Configure the logging level and optional file output for
#'   ContextNet messages.
#'
#' @param level Log level: "none", "error", "warning", "info", "debug"
#' @param file Optional file path to write logs. Set to NULL to disable
#'   file logging.
#'
#' @return The previous log level (invisibly)
#'
#' @examples
#' # Set to debug mode for troubleshooting
#' setup_logging("debug")
#'
#' # Suppress all but errors
#' setup_logging("error")
#'
#' # Log to file
#' setup_logging("info", file = "contextnet.log")
#'
#' @export
setup_logging <- function(level = c("info", "warning", "error", "debug", "none"),
                          file = NULL) {
    level <- match.arg(level)
    prev_level <- getOption("contextnet.log_level", "info")

    options(contextnet.log_level = level)

    if (!is.null(file)) {
        options(contextnet.log_file = file)
        cli::cli_alert_info("Logging to: {.path {file}}")
    }

    invisible(prev_level)
}


#' Get current log level
#'
#' @return Character string of current log level
#'
#' @keywords internal
get_log_level <- function() {
    getOption("contextnet.log_level", "info")
}


#' Internal logging function
#'
#' @param msg Message to log (supports cli formatting)
#' @param level Log level for this message
#' @param ... Additional arguments passed to cli formatting
#'
#' @return NULL invisibly
#'
#' @keywords internal
.log <- function(msg, level = "info", ...) {
    current_level <- getOption("contextnet.log_level", "info")
    levels <- c(none = 0, error = 1, warning = 2, info = 3, debug = 4)

    # Check if this message should be shown
    if (levels[level] > levels[current_level]) {
        return(invisible(NULL))
    }

    # Format the message
    formatted <- cli::format_message(msg, ...)

    # Console output
    switch(level,
        error = cli::cli_alert_danger(formatted),
        warning = cli::cli_alert_warning(formatted),
        info = cli::cli_alert_info(formatted),
        debug = cli::cli_text(cli::col_grey(paste0("[DEBUG] ", formatted)))
    )

    # File output if configured
    log_file <- getOption("contextnet.log_file")
    if (!is.null(log_file)) {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        cat(sprintf("[%s] %s: %s\n", timestamp, toupper(level), formatted),
            file = log_file, append = TRUE)
    }

    invisible(NULL)
}


#' Log a debug message
#'
#' @param msg Message to log
#' @param ... Additional arguments
#'
#' @keywords internal
.log_debug <- function(msg, ...) {
    .log(msg, level = "debug", ...)
}


#' Log an info message
#'
#' @param msg Message to log
#' @param ... Additional arguments
#'
#' @keywords internal
.log_info <- function(msg, ...) {
    .log(msg, level = "info", ...)
}


#' Log a warning message
#'
#' @param msg Message to log
#' @param ... Additional arguments
#'
#' @keywords internal
.log_warn <- function(msg, ...) {
    .log(msg, level = "warning", ...)
}


#' Log an error message
#'
#' @param msg Message to log
#' @param ... Additional arguments
#'
#' @keywords internal
.log_error <- function(msg, ...) {
    .log(msg, level = "error", ...)
}


#' Create a progress context for long operations
#'
#' @param msg Message describing the operation
#' @param verbose Logical, whether to show progress
#'
#' @return A list with start and end functions
#'
#' @keywords internal
.progress_context <- function(msg, verbose = TRUE) {
    if (!verbose) {
        return(list(
            start = function() invisible(NULL),
            end = function() invisible(NULL)
        ))
    }

    list(
        start = function() cli::cli_progress_step(msg),
        end = function() invisible(NULL)
    )
}

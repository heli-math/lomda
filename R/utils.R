# Internal helper for optional column-name arguments.
.lomda_col_arg <- function(expr, value, arg_name) {
  if (is.character(expr) && length(expr) == 1L) {
    return(expr)
  }
  if (is.symbol(expr)) {
    value <- tryCatch(value, error = function(e) NULL)
    if (is.character(value) && length(value) == 1L) {
      return(value)
    }
    return(as.character(expr))
  }
  value <- tryCatch(value, error = function(e) NULL)
  if (is.character(value) && length(value) == 1L) {
    return(value)
  }
  stop(arg_name, " must be a single column name.", call. = FALSE)
}

.lomda_bt <- function(x) {
  paste0("`", gsub("`", "\\\\`", x), "`")
}

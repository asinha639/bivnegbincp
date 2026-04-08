#' Correct Index Jumps Within a Threshold
#'
#' Cleans up a sequence of indices by merging nearby ones based on a user-defined
#' distance threshold.
#'
#' @param index A numeric vector of change point indices.
#' @param tolerance Integer. Maximum allowed difference between consecutive indices.
#'
#' @return A cleaned vector of unique indices.
#'
#' @examples
#' index_correction(c(10, 11, 12, 25, 26), tolerance = 2)
#'
#' @export
index_correction <- function(index, tolerance = 2) {
  if (!is.numeric(index)) {
    stop("Input 'index' must be a numeric vector.")
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance < 0) {
    stop("Argument 'tolerance' must be a single non-negative numeric value.")
  }

  index <- index[!is.na(index)]

  if (length(index) == 0) {
    warning("No valid index values provided. Returning empty vector.")
    return(integer(0))
  }

  index <- sort(index)

  if (length(index) > 1) {
    for (i in 2:length(index)) {
      if (abs(index[i] - index[i - 1]) <= tolerance) {
        index[i] <- index[i - 1]
      }
    }
  }

  cleaned_index <- unique(index)
  cleaned_index <- cleaned_index[cleaned_index != 0]

  if (length(cleaned_index) == 0) {
    warning("All indices removed after filtering. Returning empty vector.")
  }

  cleaned_index
}

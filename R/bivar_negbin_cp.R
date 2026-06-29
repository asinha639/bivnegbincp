#' Bivariate Negative Binomial Single Change Point Detection
#'
#' Detects a single change point in bivariate count data assumed to follow a
#' bivariate Negative Binomial model. The user chooses the prior through the
#' `prior` argument. This function serves as the core single-change-point
#' routine used by the binary segmentation procedure for multiple change point
#' detection.
#'
#' @param x1 Integer-valued numeric vector for the first count series.
#' @param x2 Integer-valued numeric vector for the second count series.
#' @param prior Prior specification. Must be either `"constant"` or
#'   `"dirichlet"`. Capitalization is ignored.
#' @param kappa_const Positive numeric value for the fixed dispersion parameter.
#' @param v1 Positive numeric hyperparameter for the Dirichlet prior.
#' @param v2 Positive numeric hyperparameter for the Dirichlet prior.
#' @param v3 Positive numeric hyperparameter for the Dirichlet prior.
#'
#' @return A one-row data frame with columns `post_prob`, `cp_index`, and
#'   `cp_status`. The value `post_prob` is the posterior probability at the
#'   selected change point, `cp_index` is the estimated change point location,
#'   and `cp_status = 1` indicates that the selected location is not too close
#'   to the boundary.
#'
#' @examples
#' set.seed(123)
#' x1 <- c(rpois(20, 8), rpois(20, 20))
#' x2 <- c(rpois(20, 2), rpois(20, 8))
#' bivnegbin_singlecp(x1, x2, prior = "constant", kappa_const = 1)
#'
#' @export
bivnegbin_singlecp <- function(x1,
                               x2,
                               prior = "constant",
                               kappa_const = 1,
                               v1 = 2,
                               v2 = 1,
                               v3 = 1) {
  prior <- .resolve_prior(prior)
  
  .validate_bnb_inputs(
    x1 = x1,
    x2 = x2,
    kappa_const = kappa_const,
    prior = prior,
    v1 = v1,
    v2 = v2,
    v3 = v3
  )
  
  if (prior == "constant") {
    return(.bivnegbin_singlecp_constant(
      x1 = x1,
      x2 = x2,
      kappa_const = kappa_const
    ))
  }
  
  .bivnegbin_singlecp_dirichlet(
    x1 = x1,
    x2 = x2,
    v1 = v1,
    v2 = v2,
    v3 = v3,
    kappa_const = kappa_const
  )
}

# Internal helper functions ------------------------------------------------

.validate_bnb_inputs <- function(x1, x2, kappa_const, prior, v1, v2, v3) {
  if (!is.numeric(x1) || !is.numeric(x2)) {
    stop("'x1' and 'x2' must be numeric vectors.")
  }
  
  if (length(x1) != length(x2)) {
    stop("Length of 'x1' and 'x2' must be equal.")
  }
  
  if (length(x1) < 4) {
    stop("At least 4 observations are required in each series.")
  }
  
  if (any(is.na(x1)) || any(is.na(x2))) {
    stop("Missing values are not allowed in 'x1' or 'x2'.")
  }
  
  if (any(x1 < 0) || any(x2 < 0)) {
    stop("All observations must be non-negative.")
  }
  
  if (any(abs(x1 - round(x1)) > .Machine$double.eps^0.5) ||
      any(abs(x2 - round(x2)) > .Machine$double.eps^0.5)) {
    stop("All observations must be integer-valued counts.")
  }
  
  if (!is.numeric(kappa_const) || length(kappa_const) != 1 || is.na(kappa_const) ||
      kappa_const <= 0) {
    stop("'kappa_const' must be a single positive number.")
  }
  
  if (prior == "dirichlet") {
    if (any(!is.numeric(c(v1, v2, v3))) ||
        any(is.na(c(v1, v2, v3))) ||
        any(c(v1, v2, v3) <= 0)) {
      stop("For prior = 'dirichlet', 'v1', 'v2', and 'v3' must be positive numbers.")
    }
  }
  
  invisible(TRUE)
}

.resolve_prior <- function(prior) {
  allowed <- c("constant", "dirichlet")
  
  if (!is.character(prior) || length(prior) != 1 || is.na(prior)) {
    stop("`prior` must be a single character string: \"constant\" or \"dirichlet\".")
  }
  
  prior_clean <- tolower(trimws(prior))
  
  if (prior_clean %in% allowed) {
    return(prior_clean)
  }
  
  d <- utils::adist(prior_clean, allowed)
  closest <- allowed[which.min(d)]
  
  if (min(d) <= 2) {
    stop(sprintf(
      "Unknown prior '%s'. Did you mean '%s'?",
      prior, closest
    ))
  }
  
  stop(sprintf(
    "Unknown prior '%s'. Please choose one of: %s.",
    prior, paste(shQuote(allowed), collapse = ", ")
  ))
}

.bivnegbin_singlecp_constant <- function(x1, x2, kappa_const = 1) {
  res <- bnb_singlecp_constant_cpp(
    x1 = as.numeric(x1),
    x2 = as.numeric(x2),
    kappa_const = kappa_const
  )

  if (is.na(res$cp_index) || is.na(res$max_post_prob)) {
    return(NA)
  }

  data.frame(
    post_prob = res$max_post_prob,
    cp_index = res$cp_index,
    cp_status = res$cp_status
  )
}

.bivnegbin_singlecp_dirichlet <- function(x1, x2, v1 = 2, v2 = 1, v3 = 1, kappa_const = 1) {
  res <- bnb_singlecp_dirichlet_cpp(
    x1 = as.numeric(x1),
    x2 = as.numeric(x2),
    v1 = v1,
    v2 = v2,
    v3 = v3,
    kappa_const = kappa_const
  )

  if (is.na(res$cp_index) || is.na(res$max_post_prob)) {
    return(NA)
  }

  data.frame(
    post_prob = res$max_post_prob,
    cp_index = res$cp_index,
    cp_status = res$cp_status
  )
}

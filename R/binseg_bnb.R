#' Binary Segmentation for Bivariate Negative Binomial Data
#'
#' Applies the bivariate Negative Binomial single change point detector
#' recursively to identify multiple change points through binary segmentation.
#' The prior is selected through the `prior` argument and the matching internal
#' single-change routine is used at each step.
#'
#' @param CGH A two-column matrix or data frame of non-negative integer counts.
#' @param prior Character string. Either `"constant"` or `"dirichlet"`.
#' @param th_cp Numeric threshold for retaining detected change points.
#' @param k_val Positive numeric fixed value of the dispersion parameter.
#' @param v1 Positive numeric hyperparameter for the Dirichlet prior.
#' @param v2 Positive numeric hyperparameter for the Dirichlet prior.
#' @param v3 Positive numeric hyperparameter for the Dirichlet prior.
#' @param save_output Logical. If `TRUE`, writes the returned table to an Excel file.
#' @param file_name Output Excel file name when `save_output = TRUE`.
#'
#' @return A data frame with columns `Posterior Probability` and `CP index`.
#' Returns `NA` when no change point exceeds the threshold.
#'
#' @examples
#' set.seed(123)
#' x1 <- c(rpois(20, 8), rpois(15, 20), rpois(15, 10))
#' x2 <- c(rpois(20, 2), rpois(15, 8), rpois(15, 3))
#' dat <- cbind(x1, x2)
#' BinSeg_BNB(dat, prior = "constant", th_cp = 0.15, k_val = 1)
#'
#' @export
BinSeg_BNB <- function(CGH,
                       prior = "constant",
                       th_cp = 0.5,
                       k_val = 5,
                       v1 = 2,
                       v2 = 1,
                       v3 = 1,
                       save_output = FALSE,
                       file_name = "bivnegbincp_output.xlsx") {
  prior <- .resolve_prior(prior)
  
  if (!is.matrix(CGH) && !is.data.frame(CGH)) {
    stop("'CGH' must be a matrix or data frame with exactly two columns.")
  }
  
  CGH <- as.data.frame(CGH)
  
  if (ncol(CGH) != 2) {
    stop("'CGH' must have exactly two columns.")
  }
  
  if (nrow(CGH) < 4) {
    stop("'CGH' must have at least 4 rows.")
  }
  
  if (!is.numeric(th_cp) || length(th_cp) != 1 || is.na(th_cp) ||
      th_cp < 0 || th_cp > 1) {
    stop("'th_cp' must be a single number between 0 and 1.")
  }
  
  index_pos_mat <- data.frame(V1 = NA_real_, V2 = NA_real_)
  new_matrix <- data.frame(V1 = numeric(0), V2 = numeric(0))
  output_frame <- data.frame(V1 = numeric(0), V2 = numeric(0))
  
  rownum <- nrow(CGH)
  level <- 0
  cpnum <- 1
  index_pos_mat[1, 1] <- 1
  index_pos_mat[1, 2] <- rownum
  
  while (cpnum > 0) {
    iter_val <- 0
    
    for (i in seq_len(cpnum)) {
      index_l <- index_pos_mat[i, 1]
      index_r <- index_pos_mat[i, 2]
      
      if ((index_r - index_l) > 2) {
        data1 <- CGH[index_l:index_r, 1]
        data2 <- CGH[index_l:index_r, 2]
        
        res <- bivnegbin_singlecp(
          x1 = data1,
          x2 = data2,
          prior = prior,
          kappa_const = k_val,
          v1 = v1,
          v2 = v2,
          v3 = v3
        )
        
        if (!anyNA(res)) {
          cp_loc <- res$cp_index
          p_value <- res$post_prob
          
          if (res$cp_status == 1) {
            if (index_l != 1) {
              cp_loc <- index_l + cp_loc - 1
            }
            
            level <- level + 1
            output_frame[level, 1] <- as.numeric(p_value)
            output_frame[level, 2] <- cp_loc
            
            iter_val <- iter_val + 1
            new_matrix[iter_val, 1] <- index_l
            new_matrix[iter_val, 2] <- cp_loc
            
            iter_val <- iter_val + 1
            new_matrix[iter_val, 1] <- cp_loc + 1
            new_matrix[iter_val, 2] <- index_r
          }
        }
      }
    }
    
    cpnum <- iter_val
    index_pos_mat <- stats::na.omit(new_matrix)
    new_matrix <- data.frame(V1 = numeric(0), V2 = numeric(0))
  }
  
  output_frame <- stats::na.omit(output_frame)
  
  if (nrow(output_frame) > 0) {
    if (max(output_frame[, 1]) < th_cp) {
      return(NA)
    }
    
    out <- output_frame[output_frame[, 1] > th_cp, , drop = FALSE]
    out <- out[order(out[, 2]), , drop = FALSE]
    names(out) <- c("Posterior Probability", "CP index")
    
    if (save_output) {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "ChangePoints")
      openxlsx::writeData(wb, sheet = "ChangePoints", x = out)
      openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
    }
    
    return(out)
  }
  
  NA
}

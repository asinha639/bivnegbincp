.reference_log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.reference_status <- function(cp_index, n) {
  if (cp_index %in% c(1L, 2L, n - 1L, n)) {
    return(0L)
  }

  1L
}

.reference_singlecp_from_log_post <- function(log_post, n) {
  post_prob <- exp(log_post - .reference_log_sum_exp(log_post))
  cp_index <- which.max(post_prob)

  data.frame(
    post_prob = max(post_prob),
    cp_index = cp_index,
    cp_status = .reference_status(cp_index, n)
  )
}

.reference_constant_segment <- function(sx1, sx2, len, kappa_const) {
  a <- len * kappa_const

  -log(a + 1) +
    lgamma(sx1 + a + 2) +
    lgamma(sx2 + 1) -
    lgamma(sx2 + sx1 + a + 3) +
    lgamma(sx1 + 1) +
    lgamma(a + 2) -
    lgamma(a + sx1 + 2)
}

.reference_dirichlet_segment_first <- function(sx1, sx2, len, v1, v2, v3, kappa_const) {
  a <- len * kappa_const

  -log(a + v3) +
    lgamma(sx1 + a + v1 + v3) +
    lgamma(sx2 + v2) -
    lgamma(sx2 + sx1 + a + v1 + v2 + v3) +
    lgamma(sx1 + v1) +
    lgamma(a + v3 + 1) -
    lgamma(a + sx1 + v1 + v3)
}

.reference_dirichlet_segment_second <- function(sx1, sx2, len, v1, v2, v3, kappa_const) {
  a <- len * kappa_const

  -log(a + v3) +
    lgamma(sx1 + a + v1 + v3) +
    lgamma(sx2 + v2) -
    lgamma(sx2 + sx1 + a + v1 + v2 + v3) +
    lgamma(sx1 + v1) +
    lgamma(a + v3 - v1 + 1) -
    lgamma(a + sx1 + v3)
}

.reference_singlecp_constant <- function(x1, x2, kappa_const = 1) {
  n <- length(x1)
  cs_x1 <- cumsum(x1)
  cs_x2 <- cumsum(x2)
  total_x1 <- cs_x1[n]
  total_x2 <- cs_x2[n]
  log_post <- numeric(n - 1)

  for (tau in seq_len(n - 1)) {
    sx11 <- cs_x1[tau]
    sx21 <- cs_x2[tau]
    sx12 <- total_x1 - sx11
    sx22 <- total_x2 - sx21

    log_post[tau] <-
      .reference_constant_segment(sx11, sx21, tau, kappa_const) +
      .reference_constant_segment(sx12, sx22, n - tau, kappa_const)
  }

  .reference_singlecp_from_log_post(log_post, n)
}

.reference_singlecp_dirichlet <- function(x1, x2, v1 = 2, v2 = 1, v3 = 1, kappa_const = 1) {
  n <- length(x1)
  cs_x1 <- cumsum(x1)
  cs_x2 <- cumsum(x2)
  total_x1 <- cs_x1[n]
  total_x2 <- cs_x2[n]
  log_post <- numeric(n - 1)

  for (tau in seq_len(n - 1)) {
    sx11 <- cs_x1[tau]
    sx21 <- cs_x2[tau]
    sx12 <- total_x1 - sx11
    sx22 <- total_x2 - sx21

    log_post[tau] <-
      .reference_dirichlet_segment_first(sx11, sx21, tau, v1, v2, v3, kappa_const) +
      .reference_dirichlet_segment_second(sx12, sx22, n - tau, v1, v2, v3, kappa_const)
  }

  .reference_singlecp_from_log_post(log_post, n)
}

.reference_binseg_bnb <- function(CGH,
                                  prior = "constant",
                                  th_cp = 0.5,
                                  k_val = 5,
                                  v1 = 2,
                                  v2 = 1,
                                  v3 = 1) {
  CGH <- as.data.frame(CGH)
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

        if (prior == "constant") {
          res <- .reference_singlecp_constant(data1, data2, kappa_const = k_val)
        } else {
          res <- .reference_singlecp_dirichlet(data1, data2, v1, v2, v3, k_val)
        }

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
    return(out)
  }

  NA
}

test_that("optimized posterior probabilities are normalized", {
  x1 <- c(8, 9, 7, 10, 22, 24, 20, 21, 11, 12)
  x2 <- c(2, 3, 2, 4, 8, 7, 9, 8, 3, 4)

  constant <- bnb_singlecp_constant_cpp(x1, x2, kappa_const = 1.25)
  dirichlet <- bnb_singlecp_dirichlet_cpp(x1, x2, v1 = 2, v2 = 1, v3 = 1, kappa_const = 1.5)

  expect_true(isTRUE(all.equal(sum(constant$post_prob), 1, tolerance = 1e-12)))
  expect_true(isTRUE(all.equal(sum(dirichlet$post_prob), 1, tolerance = 1e-12)))
  expect_equal(constant$cp_index, which.max(constant$post_prob))
  expect_equal(dirichlet$cp_index, which.max(dirichlet$post_prob))
})

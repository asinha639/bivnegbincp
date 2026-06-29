test_that("constant prior optimized backend matches reference posterior and MAP", {
  x1 <- c(8, 9, 7, 10, 22, 24, 20, 21, 11, 12)
  x2 <- c(2, 3, 2, 4, 8, 7, 9, 8, 3, 4)

  expected <- .reference_singlecp_constant(x1, x2, kappa_const = 1.25)
  observed <- bivnegbin_singlecp(x1, x2, prior = "constant", kappa_const = 1.25)

  expect_true(isTRUE(all.equal(observed$post_prob, expected$post_prob, tolerance = 1e-12)))
  expect_equal(observed$cp_index, expected$cp_index)
  expect_equal(observed$cp_status, expected$cp_status)
})

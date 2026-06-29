test_that("Dirichlet prior optimized backend matches reference posterior and MAP", {
  x1 <- c(6, 7, 8, 7, 18, 19, 17, 20, 9, 8)
  x2 <- c(3, 3, 4, 2, 9, 10, 8, 9, 4, 3)

  expected <- .reference_singlecp_dirichlet(
    x1,
    x2,
    v1 = 2,
    v2 = 1,
    v3 = 1,
    kappa_const = 1.5
  )
  observed <- bivnegbin_singlecp(
    x1,
    x2,
    prior = "dirichlet",
    v1 = 2,
    v2 = 1,
    v3 = 1,
    kappa_const = 1.5
  )

  expect_true(isTRUE(all.equal(observed$post_prob, expected$post_prob, tolerance = 1e-12)))
  expect_equal(observed$cp_index, expected$cp_index)
  expect_equal(observed$cp_status, expected$cp_status)
})

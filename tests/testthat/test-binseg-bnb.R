test_that("binary segmentation preserves change point output for constant prior", {
  x1 <- c(8, 9, 8, 10, 24, 23, 25, 21, 11, 12, 10, 13)
  x2 <- c(2, 3, 2, 4, 8, 9, 8, 7, 3, 4, 3, 5)
  dat <- cbind(x1, x2)

  expected <- .reference_binseg_bnb(dat, prior = "constant", th_cp = 0.1, k_val = 1.25)
  observed <- BinSeg_BNB(dat, prior = "constant", th_cp = 0.1, k_val = 1.25)

  expect_true(isTRUE(all.equal(observed, expected, tolerance = 1e-12)))
})

test_that("binary segmentation preserves change point output for Dirichlet prior", {
  x1 <- c(6, 7, 8, 7, 18, 19, 17, 20, 9, 8, 9, 10)
  x2 <- c(3, 3, 4, 2, 9, 10, 8, 9, 4, 3, 4, 5)
  dat <- cbind(x1, x2)

  expected <- .reference_binseg_bnb(
    dat,
    prior = "dirichlet",
    th_cp = 0.1,
    k_val = 1.5,
    v1 = 2,
    v2 = 1,
    v3 = 1
  )
  observed <- BinSeg_BNB(
    dat,
    prior = "dirichlet",
    th_cp = 0.1,
    k_val = 1.5,
    v1 = 2,
    v2 = 1,
    v3 = 1
  )

  expect_true(isTRUE(all.equal(observed, expected, tolerance = 1e-12)))
})

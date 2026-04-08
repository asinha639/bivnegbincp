test_that("bivar_negbin_cp returns a list", {
  set.seed(123)
  
  x1 <- c(
    rnbinom(20, size = 10, mu = 5),
    rnbinom(20, size = 10, mu = 12)
  )
  
  x2 <- c(
    rnbinom(20, size = 8, mu = 4),
    rnbinom(20, size = 8, mu = 10)
  )
  
  res <- bivnegbin_singlecp(
    x1 = x1,
    x2 = x2
  )
  
  expect_type(res, "list")
  expect_true(length(res) > 0)
})

test_that("bivar_negbin_cp errors when lengths differ", {
  x1 <- 1:10
  x2 <- 1:9
  
  expect_error(
    bivnegbin_singlecp(
      x1 = x1,
      x2 = x2
    )
  )
})
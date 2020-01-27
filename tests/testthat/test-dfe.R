test_that("dfgm_martin = 0 if s > so & dfgm_martin > 0 if s <= so ", {
  n <- floor(abs(rnorm(1))) + 1
  lambda <- abs(rnorm(1))
  so <- lambda * 5
  s1 <- so * 1.1
  s2 <- so * 0.9
  expect_equal(dfgm_martin(s1, n, lambda, so), 0)
  expect_gt(dfgm_martin(s2, n, lambda, so), 0)
})


test_that("dfgm_tenaillon = 0 if s > so & dfgm_tenaillon > 0 if s <= so ", {
  n      <- floor(abs(rnorm(1))) + 1
  lambda <- abs(rnorm(1))
  so     <- round(lambda * 5, 1)
  alpha  <- 1/2
  Q      <- 2
  s <- sample(seq(-2, so + 0.5, 0.1))
  fs <- dfgm_tenaillon(s, n, lambda, so, alpha, Q)
  expect_equal(fs[which(s > so)], numeric(length(which(s > so))))
  checkmate::expect_numeric(fs[which(s <= so)], lower = 0, any.missing = F, null.ok = FALSE)
})

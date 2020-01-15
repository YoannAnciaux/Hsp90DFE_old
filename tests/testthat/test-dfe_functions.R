test_that("dfgm_martin = 0 if s > so & dfgm_martin > 0 if s <= so ", {
  n <- floor(abs(rnorm(1))) + 1
  lambda <- abs(rnorm(1))
  so <- lambda * 5
  s1 <- so * 1.1
  s2 <- so * 0.9
  expect_equal(dfgm_martin(s1, n, lambda, so), 0)
  expect_gt(dfgm_martin(s2, n, lambda, so), 0)
})

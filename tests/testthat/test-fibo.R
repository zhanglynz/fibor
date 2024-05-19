test_that("fibo works", {
  v <- c(0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55)
  expect_equal(vapply(0:10, fibo, numeric(1)), v)
})

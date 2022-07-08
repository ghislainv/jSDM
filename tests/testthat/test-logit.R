#context("test-logit")

test_that("logit works", {
  expect_equal(logit(0.5), 0)
  expect_equal(logit(0.6), 0.405465108)
})

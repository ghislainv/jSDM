context("test-inv_logit")

test_that("inv_logit works", {
  expect_equal(inv_logit(0.4054651), 0.6)
  expect_equal(inv_logit(0), 0.5)
})
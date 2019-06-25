context("test-logit")

test_that("logit and inv_logit works", {
  expect_equal(logit(0.5), 0)
  expect_equal(logit(0.6), 0.4054651)
  expect_equal(inv_logit(0.4054651), 0.6)
  expect_equal(inv_logit(0), 0.5)
})

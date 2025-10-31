test_that("coalesce_num handles NULL and NA", {
  expect_equal(evolveTrial:::coalesce_num(NULL, 5), 5)
  expect_equal(evolveTrial:::coalesce_num(NA, 5), 5)
  expect_equal(evolveTrial:::coalesce_num(10, 5), 10)
  expect_equal(evolveTrial:::coalesce_num(0, 5), 0)
})

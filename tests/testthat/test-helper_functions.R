context("Helper functions")

test_that("geometric mean calculation works",{

  empty_vector <- vector()
  test_vector <- (seq(10,100))
  test_vector_with_NAs <- (rep(c(seq(20,40),NA),10))
  non_numeric_vector <- c(1,2,3,"A")
  expect_error(gm_mean(empty_vector))
  expect_error(gm_mean(non_numeric_vector))
  expect_error(gm_mean(test_vector.na.rm="TRUE"))
  expect_equal(gm_mean(5),5)
  expect_equal(round(gm_mean(test_vector)),round(47.29746))
  expect_true(is.na(gm_mean(test_vector_with_NAs,na.rm = FALSE)))
  expect_false(is.na(gm_mean(test_vector_with_NAs,na.rm = TRUE)))
})

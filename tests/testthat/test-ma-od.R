library(robustMAdesigns)

context("design matrix")

test_that("desMatrix returns a matrix with 2 rows and 5 columns for 1x2 CR design", {
  expect_equal(dim(desMatrix("1x2", "CR")), c(2, 5))
})

test_that("desMatrix returns a matrix with 2 rows and 4 columns for 1x2 CL design", {
  expect_equal(dim(desMatrix("1x2", "CL")), c(2, 4))
})

test_that("desMatrix returns a matrix with 2 rows and 4 columns for 1x2 DS design", {
  expect_equal(dim(desMatrix("1x2", "DS")), c(2, 4))
})


test_that("desMatrix returns a matrix with 12 rows and 6 columns for 2x2 DS design", {
  expect_equal(dim(desMatrix("2x2", "DS")), c(12, 6))
})


test_that("desMatrix returns a matrix with 30 rows and 8 columns for 3x2 DS design", {
  expect_equal(dim(desMatrix("3x2", "DS")), c(30, 8))
})


context("contrast matrix")

test_that("contMatrix returns a matrix with 3 rows and 5 columns for
          1x3 layout, all-pair comparison", {
  expect_equal(dim(contMatrix("1x3", "all-pair", F)), c(3, 5))
})

test_that("contMatrix returns a matrix with 3 rows and 5 columns for
          1x3 layout, global comparison", {
  expect_equal(dim(contMatrix("1x3", "global", F)), c(3, 5))
})

test_that("Construct permutation table", {

  # Note: This is a unit test to ensure `gtools::permutations()` routine
  # that we depend upon does not change.

  # K = 2, 2! x 2
  k_2 = matrix(c(1, 2, 2, 1), nrow = 2)
  expect_equal(attribute_permutation_table(2), k_2)

  # K = 3, 3! x 3
  k_3 = cbind(
    c(1, 1, 2, 2, 3, 3),
    c(2, 3, 1, 3, 1, 2),
    c(3, 2, 3, 1, 2, 1)
  )
  expect_equal(attribute_permutation_table(3), k_3)
})

test_that("Location of highest element", {
  x = matrix(1:4, nrow = 2)
  expect_equal(matrix_max_index(x), t(c(2, 2)))

  y = matrix(c(1, 2, 3, 4, 4, 5, 3, 2, 1), nrow = 3)
  expect_equal(matrix_max_index(y), t(c(3, 2)))

  z = matrix(c(1, 2, 3, 4, 5, 3, 2, 1), nrow = 2)
  expect_equal(matrix_max_index(z), t(c(1, 3)))
})

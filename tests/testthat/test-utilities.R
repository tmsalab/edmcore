## Melt Matrices and Arrays ----

test_that("Verify matrix is melted correctly", {

  x = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)

  # Construct the target data.frame
  target_df = data.frame(
    J = rep(seq_len(nrow(x)), times = ncol(x)),
    K = rep(seq_len(ncol(x)), each = nrow(x)),
    Value = c(x)
  )

  # Note: dimnames differ between the two.
  expect_equal(melt_array(x), target_df, check.attributes = FALSE)
})

test_that("Verify array is melted correctly", {

  x = array(c(1, 2, 3, 4, 5, 6), dim = c(4, 3, 2))

  # Construct the target data.frame for an array
  target_df = data.frame(
    J = rep(seq_len(nrow(x)), times = ncol(x)),
    K = rep(seq_len(ncol(x)), each = nrow(x)),
    S = rep(seq_len(dim(x)[3]), each = ncol(x) * nrow(x)),
    Value = c(x)
  )

  # Note: dimnames differ between the two.
  expect_equal(melt_array(x), target_df, check.attributes = FALSE)
})


test_that("Verify array is melted correctly", {

  x = array(c(1, 2, 3, 4, 5, 6), dim = c(4, 3, 2))

  # Construct the target data.frame for an array
  target_df = data.frame(
    J = rep(seq_len(nrow(x)), times = ncol(x)),
    K = rep(seq_len(ncol(x)), each = nrow(x)),
    S = rep(seq_len(dim(x)[3]), each = ncol(x) * nrow(x)),
    Value = c(x)
  )

  # Note: dimnames differ between the two.
  expect_equal(melt_array(x), target_df, check.attributes = FALSE)
})

## List Matrices ----

test_that("Verify list of matrices coerces to array", {

  # Build a list
  x = matrix(1:6, nrow = 2)
  x_list = list(x, x + 1)

  # Build expected array
  target_array = array(c(x, x+1), dim = c(2, 3, 2))

  # Verify equality
  expect_equal(listmatrix_to_array(x_list), target_array)
})

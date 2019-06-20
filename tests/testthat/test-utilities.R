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


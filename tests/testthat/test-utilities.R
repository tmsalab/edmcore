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

## Sequence Generation ----

test_that("Sequence generation decreasing", {

  # Generate sequence for k = 3
  k = 3

  # (k-1):-1:0 => 2, 1, 0
  # Details:https://www.mathworks.com/help/matlab/ref/colon.html
  r_series = t(seq(k - 1, 0, by = -1))
  cpp_series = seq_linear_decrease(k - 1, 0)

  expect_equal(r_series, cpp_series)
})

test_that("Sequence generation increasing", {

  # Generate sequence for k = 3
  k = 3

  # (k-1):-1:0 => 2, 1, 0
  # Details:https://www.mathworks.com/help/matlab/ref/colon.html
  r_series = t(seq(0, k - 1, by = 1))
  cpp_series = seq_linear_increase(0, k - 1)

  expect_equal(r_series, cpp_series)
})

## Binary Ideal Q Matrix ----

test_that("Compute the binary matrix", {

  r_binary_matrix = function(k) {
    # Required
    if(k == 1)
      return(matrix(1))

    # Construct matrix with identical columns
    x = replicate(k, 1:(2^k - 1))

    # Construct a row vector
    division_b = 2^((k-1):0)

    # Element-wise division of matrix rows.
    divs =  t(t(x) / division_b)

    # Floor the results
    divs = floor(divs)

    divs[, 2:k] = divs[, 2:k] - 2*divs[, 1:(k-1)]

    divs
  }

  # Test special case
  # Generate sequence for k = 1
  k = 1
  expect_equal(r_binary_matrix(k), binary_q_ideal(k))

  # Generate sequence for k = 2
  k = 2
  expect_equal(r_binary_matrix(k), binary_q_ideal(k))

  # Generate sequence for k = 3
  k = 3
  expect_equal(r_binary_matrix(k), binary_q_ideal(k))

  # Generate sequence for k = 2
  k = 4
  expect_equal(r_binary_matrix(k), binary_q_ideal(k))
})

## Combinations ----
test_that("Generate All Combinations", {
  skip_if_not_installed("gtools")

  expect_equal(edmcore:::combination_matrix(5, 2),
               gtools::combinations(5, 2))

  expect_equal(edmcore:::combination_matrix(8, 3),
               gtools::combinations(8, 3))
})

test_that("Generate Combinations under Set", {
  skip_if_not_installed("gtools")

  expect_equal(edmcore:::combination_matrix_from_vector(3:5, 2),
               gtools::combinations(3, 2, v = 3:5))

  expect_equal(edmcore:::combination_matrix_from_vector(4:7, 3),
               gtools::combinations(4, 3, v = 4:7))
})

test_that("Compute n Choose k ", {
  expect_equal(edmcore:::n_choose_k(5, 2),
               base::choose(5, 2))

  expect_equal(edmcore:::n_choose_k(8, 3),
               base::choose(8, 3))
})

## Set Difference ----

test_that("Finds missing elements", {
  x = c(1, 2, 3)
  y = c(2, 3, 4)
  expect_equal(edmcore:::set_diff_cpp(x, y),
               base::setdiff(x, y), check.attributes = FALSE)

  x = c(1, 2, 3)
  y = c(1, 4, 3, 5)
  expect_equal(edmcore:::set_diff_cpp(x, y),
               base::setdiff(x, y), check.attributes = FALSE)
})

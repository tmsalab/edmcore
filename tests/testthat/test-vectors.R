context("vectors")


# Same number of attributes
K = 5

test_that("Bijection", {

  check_val = matrix(2^((K-1):0), nrow = K, ncol = 1)

  expect_equal(bijectionvector(K),  check_val)
})


test_that("Inverse Bijection", {

  alpha = logical(K)
  cl = 5
  cl_r = cl

  for(k in (seq_len(K) - 1)) {
    pow2 = 2^(K - k - 1)
    alpha[k+1] = (pow2 <= cl_r)
    cl_r = cl_r - pow2 * alpha[k+1]
  }

  expect_equal(as.logical(inv_bijectionvector(K, cl)), alpha)
})

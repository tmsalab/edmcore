test_that("Valid attribute_bijection", {
  bijection_test = matrix( 2^((3-1):0) )

  expect_equal(attribute_bijection(3),  bijection_test)
})

test_that("Valid attribute_inv_bijection", {

  inv_bijection_test = function(K, CL) {
    alpha = rep(NA, K)

    for (k in seq_len(K) - 1) {
      twopow = 2^(K - k - 1)
      alpha[k + 1] = 1 * (twopow <= CL)
      CL = CL - twopow * alpha[k + 1]
    }

    matrix(alpha)
  }

  expect_equal(attribute_inv_bijection(3, 0),  inv_bijection_test(3, 0))

  expect_equal(attribute_inv_bijection(4, 1),  inv_bijection_test(4, 1))

  expect_equal(attribute_inv_bijection(5, 2),  inv_bijection_test(5, 2))

})

test_that("Valid attribute inverse general bijection", {
  inv_gen_bijection_test = function(K, M, CL) {
    alpha = rep(NA, K)
    for (k in (seq_len(K)-1)) {
      Mpow = M^(K - k - 1)
      ak = 0.;
      while (((ak + 1) * Mpow <= CL) & (ak < M)) {
        ak = ak + 1
      }
      alpha[k+1] = ak
      CL = CL - Mpow * alpha[k+1]
    }
    matrix(alpha)
  }

  expect_equal(attribute_inv_gen_bijection(3, 2, 0),  inv_gen_bijection_test(3, 2, 0))
  expect_equal(attribute_inv_gen_bijection(4, 2, 1),  inv_gen_bijection_test(4, 2, 1))
  expect_equal(attribute_inv_gen_bijection(5, 2, 2),  inv_gen_bijection_test(5, 2, 2))
})

test_that("Valid attribute classes (alpha matrix / pi references)", {

  # Old, r-specific pi mapping function
  pi_reference = function(K) {
    biject.vector = attribute_bijection(K)
    As = as.matrix(
      expand.grid( rep( list(c(0, 1)), K) )
    )
    a = As %*% biject.vector
    As = As[a + 1,]
    return(As)
  }

  # Check equality
  expect_equal(attribute_classes(2), pi_reference(2), check.attributes = FALSE,
               info = "Verify latent class mapping is correct.")

  # Check equality
  expect_equal(attribute_classes(5), pi_reference(5), check.attributes = FALSE,
               info = "Verify latent class mapping is correct.")

  # Check equality
  expect_equal(attribute_classes(8), pi_reference(8), check.attributes = FALSE,
               info = "Verify latent class mapping is correct.")
})

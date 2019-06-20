test_that("Verify Q matrix construction", {

  # Prototype 7 x 2
  x = matrix(c(0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0), ncol = 2)
  q = q_matrix(x)

  # Ensure class assignment
  expect_s3_class(q, c("q_matrix", "matrix"))

  # Check spacing
  expect_equal(dimnames(q),
               list(paste0("Item", seq_len(nrow(x))),
                    paste0("Trait", seq_len(ncol(x)))
               )
  )

  # Verify identifiability set
  expect_true(attr(q, "identifiable"))

  # Reset q and focus on element entries
  class(q) = "matrix"
  expect_equal(q, x, check.attributes = FALSE)
})

test_that("Verify non-identifiable Q matrix construction", {

  x = matrix(c(0, 1), ncol = 2)
  q = q_matrix(x)

  # Verify matrix isn't identifiabile
  expect_false(attr(q, "identifiable"))

  # Reset q and focus on element entries
  class(q) = "matrix"
  expect_equal(q, x, check.attributes = FALSE)
})

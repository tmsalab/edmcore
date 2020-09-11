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

  # Verify matrix is identifiable
  expect_true(is_q_strict(q))

  # Reset q and focus on element entries
  class(q) = "matrix"
  expect_equal(q, x, check.attributes = FALSE)
})

test_that("Verify non-identifiable Q matrix construction", {

  x = matrix(c(0, 1), ncol = 2)
  q = q_matrix(x)

  # Verify matrix isn't identifiable
  expect_false(is_q_strict(q))

  # Reset q and focus on element entries
  class(q) = "matrix"
  expect_equal(q, x, check.attributes = FALSE)
})


## Generic Q Matrix ----

test_that("Q is generic complete", {

  # Test a q3 matrix
  q3_generic = rbind(diag(3),
                     c(1, 1, 0),
                     c(1, 0, 1),
                     c(0, 1, 1),
                     c(1, 1, 1))

  expect_true(is_q_generic_complete(q3_generic))

  # Test a q4 matrix
  q4_generic_complete = matrix(1, nrow = 7, ncol = 4)

  expect_true(is_q_generic_complete(q4_generic_complete))

  # Try with elements zero'd out
  q3_fail = t(c(0, 0, 0))

  expect_false(is_q_generic_complete(q3_fail))

})


test_that("Q is generic identified", {
  q3_generic_id = rbind(c(0, 0, 1),
            c(1, 0, 1),
            c(1, 1, 0),
            c(0, 1, 1),
            c(1, 0, 1),
            c(1, 1, 0),
            c(1, 1, 1))

  expect_true(edmcore:::is_q_generic_identified(q3_generic_id))

  q4_generic_id = rbind(
    c(1, 1, 1, 0),
    c(1, 1, 0, 0),
    c(0, 0, 1, 1),
    c(0, 1, 0, 1),
    c(1, 0, 0, 1),
    c(0, 1, 1, 0),
    c(0, 1, 1, 1),
    c(1, 0, 0, 1),
    c(1, 1, 0, 0),
    c(0, 0, 1, 1)
  )

  expect_true(edmcore:::is_q_generic_identified(q4_generic_id))
})

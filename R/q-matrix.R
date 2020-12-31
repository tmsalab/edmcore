## Formatting Tools for Q Matrices ----

#' Format Q Matrix
#'
#' Applies a common naming scheme to the Q Matrix.
#'
#' @param x A [`base::matrix()`] with dimensions \eqn{J \times K}{J x K}.
#'
#' @return
#' A `matrix` with:
#'
#' - **Columns** named as `TraitXYZ` with XYZ denoting the trait number.
#' - **Rows** named as `ItemXYZ` with XYZ denoting the item number.
#'
#' @noRd
format_q_matrix = function(x) {

  # Extract dimensions
  j = nrow(x)
  k = ncol(x)

  # Pad naming with 0's
  colnames(x) = sprintf(paste0("Trait%0", nchar(k), "d"), seq_len(k))
  rownames(x) = sprintf(paste0("Item%0", nchar(j), "d"), seq_len(j))

  # Release
  x
}

## Verify Q matrix object ----

#' Is object a Q Matrix?
#'
#' @param x An object to test
#'
#' @return
#' A logical vector
#'
#' @export
#' @examples
#' x = matrix(c(1, 0, 1, 1))
#' is_q_matrix(q_matrix(x))
is_q_matrix = function(x) {
  inherits(x, "q_matrix")
}

#' Is Q Matrix Strictly or Generically Identifiable?
#'
#' @param x A [q_matrix()] or [base::matrix()] to test.
#'
#' @return
#' A logical vector
#'
#' @section Strictly Identifiable:
#'
#' If \eqn{\mathbf{Q}} is in the strictly identifiable set \eqn{\mathcal{Q}},
#' then it must satisfy the following conditions:
#'
#' - **(C1)** The rows of \eqn{\boldsymbol{Q}} can be permuted to the form,
#'    \eqn{\boldsymbol{Q}^\top=\left[{\boldsymbol{I_K},\boldsymbol{I_K}, (\boldsymbol{Q}^\ast)^\top}\right]^\top}
#'    where \eqn{\boldsymbol{I_K}} is a \eqn{K}-dimensional identity matrix and
#'    \eqn{\boldsymbol{Q}^\ast} is a \eqn{(J-2K)\times K} matrix.
#' - **(C2)** For any two latent classes \eqn{c} and \eqn{c'}, there exists at least one
#'    item in \eqn{\boldsymbol{Q}^\ast}, in which
#'    \eqn{\boldsymbol{\theta}_{jc}\neq \boldsymbol{\theta}_{jc'}}.
#'
#' In a more practical light, this means **(C1)** requires \eqn{\boldsymbol{Q}}
#' to include two simple structure items for each attribute and **(C2)**
#' states there must be at least one item not specified for **(C1)**
#' that distinguishes between all pairs of classes.
#'
#' @rdname is-q-identified
#' @export
#' @examples
#' x = matrix(c(1, 0, 1, 1))
#' is_q_strict(x)
#' is_q_strict(q_matrix(x))
is_q_strict = function(x) {
  stopifnot("`x` must be either `q_matrix` or `matrix`" =
              inherits(x, c("q_matrix", "matrix")))

  if (is_q_matrix(x)) {
    attr(x, 'strictly_identifiable')
  } else {
    stopifnot("`x` must contain only 0 or 1 entries." = x %in% c(0, 1))
    is_strict_q_identified(x)
  }
}

#' @section Generically Identifiable:
#'
#' Add details...
#'
#' @rdname is-q-identified
#' @export
#' @examples
#' is_q_generic(x)
#' is_q_generic(q_matrix(x))
is_q_generic = function(x) {
  stopifnot("`x` must be either `q_matrix` or `matrix`" =
              inherits(x, c("q_matrix", "matrix")))

  if (is_q_matrix(x)) {
    attr(x, 'generically_identifiable')
  } else {
    stopifnot("`x` must contain only 0 or 1 entries." = x %in% c(0, 1))
    is_q_generic_identified(x)
  }
}

## Convert Data to Q Matrix Object ----

#' Constructor for Q Matrix
#'
#' Standardizes the initialization for a Q Matrix in _R_.
#'
#' @param x A [`base::matrix()`] with dimensions \eqn{J x K}.
#'
#' @return
#' A `q_matrix` object with a fallback to `matrix`.
#'
#' @noRd
create_q_matrix = function(x) {

  # Verify strict identifiability
  strictly_identified_q = is_q_strict(x)

  generically_identified_q = if(strictly_identified_q) {
    TRUE
  } else {
    is_q_generic_identified(x)
  }

  # Structure Q matrix
  x = format_q_matrix(x)

  # Apply classes
  class(x) = c('q_matrix', class(x))

  # Embed information
  attr(x, 'strictly_identifiable') = strictly_identified_q
  attr(x, 'generically_identifiable') = generically_identified_q

  # Release result
  x
}

#' Create a Q Matrix Object
#'
#' Provides a way to create an object as a `"q_matrix"`.
#'
#' @param x        Either a `data.frame` or `matrix`.
#'
#' @return
#' A `q_matrix` object.
#'
#' @seealso
#' [as_q_matrix()]
#'
#' @export
#' @examples
#' # Q matrix values
#' x = matrix(c(1, 0, 0, 1), nrow = 2)
#'
#' # Q matrix wrapper
#' q_mat = q_matrix(x)
#'
#' # View the converted Q matrix
#' q_mat
#'
#' # Data Frame encoding of Q
#' q_df = data.frame(
#'    k1 = c(1, 0),
#'    k2 = c(0, 1)
#' )
#'
#' # Create a Q matrix
#' q_mat = q_matrix(q_df)
q_matrix = function(x) {
  as_q_matrix(x)
}

#' Coerce `data.frame` and `matrix` classes to Q Matrix.
#'
#' `as.q_matrix` acts as an aliases.
#'
#' @param x        Either a `data.frame` or `matrix`.
#' @param ...      Not used
#'
#' @return
#' A `q_matrix` object.`
#'
#' @seealso
#' [q_matrix()]
#'
#' @export
#' @rdname as_q_matrix
#'
#' @examples
#' # Q matrix values
#' x = matrix(c(1, 0, 0, 1), nrow = 2)
#'
#' # Construct class
#' q_mat = as_q_matrix(x)
#'
#' # View the converted Q matrix
#' q_mat
as_q_matrix = function(x, ...) {
  UseMethod("as_q_matrix")
}

#' @export
#' @rdname as_q_matrix
as_q_matrix.data.frame = function(x, ...) {
  x = data.matrix(x)
  create_q_matrix(x)
}

#' @export
#' @rdname as_q_matrix
as_q_matrix.matrix = function(x, ...) {
  create_q_matrix(x)
}

#' @export
#' @rdname as_q_matrix
as_q_matrix.default = function(x, ...) {
  stop(class(x)[1], " is not yet supported for conversion to `q_matrix`.")
}

## Output Q Matrix ----

#' Output Q Matrix
#'
#' Custom print method for the Q Matrix Object.
#'
#' @param x        An `q_matrix` object
#' @param ...      Additional methods passed onto the `print.matrix` method.
#'
#' @seealso
#' [q_matrix()], [as_q_matrix()]
#'
#' @return
#' An invisible `matrix` without the `q_matrix` class displayed as a part
#' of the output displayed.
#'
#' @export
#' @examples
#' # Q matrix values
#' x = matrix(c(1, 0, 0, 1), nrow = 2)
#'
#' # Show Q matrix structure
#' q_matrix(x)
print.q_matrix = function(x, ... ) {

  cat("A q matrix:", nrow(x), "x", ncol(x), "\n")
  cat("Strictly Identified: ")

  # Creative use of STDERROR vs. STDOUT. Might back fire.
  if(is_q_strict(x)) {
    cat("Yes. \n\n")
  } else {
    message("No.\n")
  }

  cat("Generically Identified: ")

  # Creative use of STDERROR vs. STDOUT. Might back fire.
  if(is_q_generic(x)) {
    cat("Yes. \n\n")
  } else {
    message("No.\n")
  }

  y = x
  class(y) = c("matrix")
  attributes(y)["strictly_identifiable"] = NULL
  attributes(y)["generically_identifiable"] = NULL
  print(y, ...)
  invisible(x)
}


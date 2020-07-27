## Formatting Tools for Item Matrices ----

#' Format Item Matrix
#'
#' Applies a common naming scheme to the Q Matrix.
#'
#' @param x A [`base::matrix()`] with dimensions \eqn{J \times K}{J x K}.
#'
#' @return
#' A `matrix` with:
#'
#' - **Rows** named as `SubjectXYZ` with XYZ denoting the subject number.
#' - **Columns** named as `ItemXYZ` with XYZ denoting the item number.
#'
#' @noRd
format_item_matrix = function(x) {

  # Extract dimensions
  n = nrow(x)
  j = ncol(x)

  # Pad naming with 0's
  rownames(x) = sprintf(paste0("Subject%0", nchar(n), "d"), seq_len(n))
  colnames(x) = sprintf(paste0("Item%0", nchar(j), "d"), seq_len(j))

  # Release
  x
}

## Verify Item Matrix Object ----

#' Is object an Item Matrix?
#'
#' @param x An object to test
#'
#' @return
#' A logical vector
#'
#' @export
#' @examples
#' # Specify int matrix details
#' x = matrix(c(1, 0, 1, 1, 1), ncol = 2)
#'
#' # Coerce to an item matrix
#' items = item_matrix(x)
#'
#' # Verify output is an item matrix
#' is_item_matrix(items)
#'
#' # Not an item matrix, but a regular matrix...
#' is_item_matrix(x)
is_item_matrix = function(x) {
  inherits(x, "item_matrix")
}

## Convert Data to an Item Matrix Object ----

#' Constructor for Item Matrix
#'
#' Standardizes the initialization for an Item Matrix in _R_.
#'
#' @param x A [`base::matrix()`] with dimensions \eqn{N \times J}{N x J}.
#'
#' @return
#' A `item_matrix` object with a fallback to `matrix`.
#'
#' @noRd
create_item_matrix = function(x) {

  # Structure Q matrix
  x = format_item_matrix(x)

  # Apply classes
  class(x) = c('item_matrix', class(x))

  # Release result
  x
}

#' Create an Item Matrix Object
#'
#' Provides a way to create an object as a `"item_matrix"`.
#'
#' @param x        Either a `data.frame` or `matrix`.
#'
#' @return
#' A `item_matrix` object.
#'
#' @seealso
#' [as_item_matrix()]
#'
#' @export
#' @examples
#' # Item matrix values
#' x = matrix(c(1, 0, 0, 1, 1, 1), ncol = 2)
#'
#' # Item matrix initialization
#' items = item_matrix(x)
#'
#' # View item matrix
#' items
#'
#' # Data Frame encoding of items
#' items_df = data.frame(
#'    s1 = c(1, 0, 1, 1, 0, 1),
#'    s2 = c(0, 1, 0, 1, 1, 0)
#' )
#'
#' # Create an Item Matrix
#' items = item_matrix(items_df)
#'
#' # View item matrix
#' items
item_matrix = function(x) {
  as_item_matrix(x)
}

#' Coerce `data.frame` and `matrix` classes to an Item Matrix.
#'
#' Converts underlying data into an Item Matrix frame.
#'
#' @param x        Either a `data.frame` or `matrix`.
#' @param ...      Not used
#'
#' @return
#' An `item_matrix` object.`
#'
#' @seealso
#' [item_matrix()]
#'
#' @export
#'
#' @details
#' [`item_matrix()`] acts as an alias for [`as_item_matrix()`]
#'
#' @examples
#' # Item matrix values
#' x = matrix(c(1, 0, 0, 1), nrow = 2)
#'
#' # Construct item class
#' items = as_item_matrix(x)
#'
#' # View output
#' items
as_item_matrix = function(x, ...) {
  UseMethod("as_item_matrix")
}

#' @export
as_item_matrix.data.frame = function(x, ...) {
  x = data.matrix(x)
  create_item_matrix(x)
}

#' @export
as_item_matrix.matrix = function(x, ...) {
  create_item_matrix(x)
}

#' @export
as_item_matrix.default = function(x, ...) {
  stop(class(x)[1], " is not yet supported for conversion to `item_matrix`.")
}

## Output Item Matrix ----

#' Output Item Matrix
#'
#' Custom print method for the Item Matrix Object.
#'
#' @param x        An `item_matrix` object
#' @param ...      Additional methods passed onto the `print.matrix` method.
#'
#' @seealso
#' [item_matrix()], [as_item_matrix()]
#'
#' @return
#' An invisible `matrix` without the `item_matrix` class displayed as a part
#' of the output displayed.
#'
#' @export
#' @examples
#' # Item matrix values
#' x = matrix(c(1, 0, 0, 1), nrow = 2)
#'
#' # Show item matrix structure
#' item_matrix(x)
print.item_matrix = function(x, ... ) {

  cat("An item matrix:", nrow(x), "x", ncol(x), "\n\n")

  # Call matrix print method
  old_class = class(x)
  class(x) = c("matrix")
  print(x, ...)
  class(x) = old_class

  invisible(x)
}


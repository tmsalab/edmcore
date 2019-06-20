#' Convert an Array into a Long data.frame
#'
#' Takes either an `array` or `matrix` and converts it into a long-form
#' `data.frame` for use with `ggplot2`.
#'
#' @param x A `matrix` or `array`.
#'
#' @return
#'
#' A `data.frame` with either 3 or 4 columns depending on whether
#' `x` is a `matrix` or `array`.
#'
#' - `J`: Row Number
#' - `K`: Column Number
#' - `S`: Depth (if `array`)
#' - `Value`: Entry in `matrix` or `array` at position.`
#'
#' @details
#' The implementation is similar to [as.data.frame.table()].
#' Though, this is customized with element naming and
#' allows us to avoid an external dependency on `reshape2`.
#'
#' @examples
#' x = matrix(1:6, nrow = 2)
#' melt_array(x)
#'
#' @noRd
melt_array = function(x) {

  stopifnot(length(dim(x)) %in% c(2L, 3L))

  j = seq_len(nrow(x))
  k = seq_len(ncol(x))

  df_matrix_long =
    if(length(dim(x)) == 2L) expand.grid(J = j,  K = k) else expand.grid(J = j,  K = k, S = seq_len(dim(x)[3]))

  # Flatten the array to a vector
  df_matrix_long$Value = c(x)

  # Return melted value
  df_matrix_long
}

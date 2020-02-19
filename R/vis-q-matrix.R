#' Heatmap Visualization of a Q Matrix
#'
#' Provides a heatmap approach to showing the estimated binary or averaged
#' values of the Q Matrix.
#'
#' @param x   Either an `edina`, `errum`, or `q_matrix` object.
#' @param ... Additional paramters not used
#'
#' @return
#' A `ggplot2` object that can be further manipulated.
#'
#' @importFrom ggplot2 ggplot aes labs theme element_text theme_minimal scale_fill_gradient geom_tile scale_x_continuous scale_y_continuous
#' @export
#' @rdname q_graph
#' @examples
#' q = q_matrix(matrix(c(1, 0, 1, 1, 0, 1), ncol = 3))
#' plot(q)
plot.q_matrix = function(x, ...) {

  # TODO: Figure out the best way to deal with globals
  Trait = Item = Value = NULL

  j = nrow(x)
  k = ncol(x)

  df_matrix_long = melt_array(x)
  colnames(df_matrix_long) = c("Item", "Trait", "Value")

  ggplot(df_matrix_long, aes(x = Trait, y = Item)) +
    geom_tile(aes(fill = Value), color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Heatmap of Q Matrix Entries",
         subtitle = sprintf("J = %i, K = %i", j, k),
         x = "Traits",
         y = "Items",
         fill = "Value") +
    theme_minimal() +
    theme() +
    scale_x_continuous(breaks = seq_len(k)) +
    scale_y_continuous(breaks = seq_len(j))
}

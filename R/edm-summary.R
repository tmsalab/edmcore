#' Construct an EDM Summary class object
#'
#' Base class that provides slots for growth within the modeling objects.
#'
#' @param property  Information on the data and algorithm configuration.
#' @param estimate  Mean, SD, Quantiles on chains pased in.
#' @param ...       Additional model-specific data
#' @param class     Specify a model-specific summary class.
#'
#' @return
#' An EDM summary object with subclasses of `class`, `edm_summary`, and `list` that
#' contains:
#'
#' - `property`: Estimation details including:
#'    - `n`: Number of Observations
#'    - `k`: Number of Attributes
#'    - `order`: Highest degree from 1 to P = 2^K or less.
#' - `estimate`: Mean, SD, Quantiles on MCMC chains.
#'
#' @export
new_edm_summary = function(property, estimate, ..., class = character()) {
  structure(
    list(
      "property" = property,
      "estimate" = estimate,
      ...
    ),
    class = c(class, "edm_summary", "list")
  )
}

## Summary S3 Extension method ----

#' Summarize MCMC chains under an EDM model
#'
#' @inherit base::summary
#'
#' @export
summary.edm_model = function(object, ...) {

  # Obtain the MCMC estimates
  mcmc_estimates = lapply(object$chain, mcmc_chain_summary_dispatch)

  # Check if we need to elevate any estimates
  if(!is.null(object$estimate)) {
    mcmc_estimates = c(mcmc_estimates, object$estimate)
  }

  new_edm_summary(
    property = object$property,
    estimate = mcmc_estimates,
    ...
  )
}


## Internal functions that compute summary statistics ----

# Obtain summary statistics on each component of the MCMC chain
summarize_mcmc = function(x, margin,
                          probs = NULL,
                          transpose_quantiles = FALSE) {
  if (!is.null(probs)) {
    quantiles = apply(
      x,
      MARGIN = margin,
      FUN = stats::quantile,
      probs = probs,
      na.rm = TRUE
    )

    if (transpose_quantiles) {
      quantiles = t(quantiles)
    }
  } else {
    quantiles = NULL
  }

  list(
    mean      = apply(x, MARGIN = margin, FUN = mean,
                      na.rm = TRUE),
    sd        = apply(x, MARGIN = margin, FUN = stats::sd,
                      na.rm = TRUE),
    quantiles = quantiles
  )
}

# Obtain an output related to the dimension of how the MCMC chain was stored.
summarize_1d_vector_output = function(x, probs) {
  summarize_mcmc(x, margin = 2,
                 probs = probs,
                 transpose_quantiles = TRUE)
}

summarize_2d_matrix_output = function(x, probs) {
  summarize_mcmc(x, margin = 1,
                 probs = probs,
                 transpose_quantiles = TRUE)
}

summarize_3d_array_output = function(x, probs) {
  summarize_mcmc(x, margin = c(1, 2),
                 probs = probs)
}

summarize_4d_array_output = function(x) {
  summarize_mcmc(x, margin = c(1, 2, 3))
}

# Route to the correct summary manipulation
mcmc_chain_summary_dispatch = function(x,
                                       probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) {
  # Retrieve properties of object
  class_x = class(x)
  type_x = typeof(x)
  dim_x = dim(x)
  n_cols = dim_x[2]
  n_dim_x = length(dim(x))

  if (class_x == "matrix" && type_x == "double" && n_cols == 1L) {
    # Handles the 1D matrix case
    summarize_1d_vector_output(x, probs = probs)
  } else if (class_x == "matrix" && type_x == "double" && n_dim_x == 2L) {
    # Handles the 2D matrix case
    summarize_2d_matrix_output(x, probs = probs)
  } else if(class_x == "array" && type_x == "double" && n_dim_x == 3L) {
    # Handles the 3D cube case
    summarize_3d_array_output(x, probs = probs)
  } else if(class_x == "array" && type_x == "double" && n_dim_x == 4L) {
    # Handles the 4D cube case
    summarize_4d_array_output(x)
  } else {
    stop("Dimensions: d = ", dim_x, " is not supported. Please try with 4 or less.")
  }
}

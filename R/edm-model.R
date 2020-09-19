#' Construct a property list for the model.
#'
#' Holds detail information about the model.
#'
#' @param n       Number of Observations
#' @param j       Number of Items
#' @param k       Number of Attributes
#' @param runtime Duration of C++ Estimation
#' @param ...     Additional properties
#' @param class   Designate the model-specific class.
#'
#' @export
new_edm_default_property_list = function(n, j, k, runtime, ..., class = character()) {
  structure(
    list(
      n = n,
      k = k,
      runtime = runtime,
      ...,
      call = match.call()
    ),
    class = c(class, "edm_property_list", "list")
  )
}

#' Display Exploratory Diagnostic Model Properties
#'
#' Custom print method for viewing Exploratory Diagnostic Model properties.
#' @inherit base::print
#' @export
print.edm_property_list = function(x, ...) {
  # Display call
  cat("Call:\n")
  print(x$call)

  # Supress call
  public_features = x
  public_features$call = NULL

  # Show the model level information
  cat("Model details:\n")
  cat(paste0("* ", names(x), ": ", x, collapse = "\n"))

  invisible(x)
}

#' Construct an EDM model class object with pre-defined names.
#'
#' Base class that provides slots for growth within the modeling objects.
#'
#' @param chain     Values from the MCMC chains across iterations.
#' @param property  Information on the data and algorithm configuration.
#' @param estimate  Any estimate done within C++ and raised into R during the
#'                  modeling procedure.
#' @param ...       Additional model-specific data
#'
#' @return
#' An EDM model object with subclasses of `class`, `edm`, and `list` that
#' contains:
#'
#' - `chains`: Chain values
#' - `property`: Estimation details including:
#'    - `n`: Number of Observations
#'    - `k`: Number of Attributes
#'    - `order`: Highest degree from 1 to P = 2^K or less.
#' - `estimates`:
#'    - **Note:** May be `NA` or a sub-list.
#'
#' @export
new_edm_model = function(chain, property_list, estimate = NA, ..., class = character()) {
  structure(
    list(
      "chain" = chain,
      "property_list" = property_list,
      "estimate" = estimate,
      ...
    ),
    class = c(class, "edm_model", "list")
  )
}

#' Print Exploratory Diagnostic Model Components
#'
#' Custom print method for viewing Exploratory Diagnostic Model components.
#' @inherit base::print
#' @export
print.edm_model = function(x, ...) {
  print(x$property_list)
  print(x$estimate)
}

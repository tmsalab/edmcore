#' Create a Generic Matrix Read Function
#'
#' Imports data from a flat file as a matrix with a predefined
#' naming scheme.
#'
#' @return Provides a function that can be customized.
#' @keywords internal
#' @importFrom utils read.table
#'
#' @noRd
read_psych = function(file, header = FALSE, sep = " ", skip = 0) {

  # Read in data as a data.frame
  read.table(file, header = header, sep = sep, skip = skip)

}


#' Import an Items Matrix
#'
#' Allows for a dichotomous items matrix to be imported with standard styling.
#'
#' @param file         name, url, or [connections] of the file to be read
#' @param header       a logical value to indicate if the data contains a
#'                     description field as the first row. Default is `FALSE`.
#' @param sep          the field separator character. Values on each line of the
#'                     file are separated by this character. Default is white
#'                     space given by `sep = " "`.
#' @param skip         the number of lines of the data file that should be
#'                     skipped before the data is read. Default is 0.
#'
#' @return
#' A `matrix` labeled with row names as `SubjectXX` and column names as `ItemYY`
#'
#' @seealso
#' [item_matrix()] / [as_item_matrix()] for constructing an Item Matrix from an
#' existing class.
#'
#' @export
#' @details
#' This function is designed to specifically read in dichotomous item matrices
#' into _R_. The matrix must be structured with the appropriate separator.
read_item_matrix = function(file, header = FALSE, sep = " ", skip = 0) {
  as_item_matrix(read_psych(file, header = header, sep = sep, skip = skip))
}

#' Import a Q Matrix
#'
#' Allows for a dichotomous Q Matrix to be imported with standard styling.
#'
#' @inheritParams read_item_matrix
#'
#' @return
#' A `matrix` labeled with row names as `ItemYY` and column names as `SkillZZ`
#'
#' @seealso
#' [q_matrix()] or [as_q_matrix()] for constructing a Q matrix from an
#' existing class.
#'
#' @export
#' @details
#' This function is designed to specifically read in dichotomous Q matrix
#' into _R_. The matrix must be structured with the appropriate separator.
read_q_matrix = function(file, header = FALSE, sep = " ", skip = 0) {
  as_q_matrix(read_psych(file, header = header, sep = sep, skip = skip))
}

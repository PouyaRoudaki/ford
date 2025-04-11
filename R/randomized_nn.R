#' Shuffle a Vector If It Has More Than One Element
#'
#' This helper function shuffles a vector only if it contains more than one element.
#' If the input is a single value (length 1), it simply returns it unchanged.
#'
#' @param x A vector to shuffle.
#'
#' @return A shuffled version of the input vector, or the original value if `x` has length 1.
#' @keywords internal
#'
#' @examples
#' .shuffle(c(1, 2, 3))   # Might return c(3, 1, 2)
#' .shuffle(7)            # Returns 7
.shuffle <- function(x) {
  if (length(x) > 1) sample(x) else x
}


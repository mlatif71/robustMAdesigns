#' @title Estimability of an effect
#'
#' @description Check the estimability of an effect
#'
#' @return A logical variable indicating whether the given contrast is
#' estimable
#'
#' @export estimable
#'
#' @param dmat design matrix of interest
#'
#' @param cmat contrast matrix corresponding to the effect of interest
#'
#' @references Searle SR (1971). Linear Models. Wiley.
#'
#' @seealso \code{contMatrix}, \code{desMatrix}
#'
#' @examples
#' # design matrix
#'   dmat <- desMatrix("1x3", design = "CR")
#'
#' # contraxt vector for comparing first and second factor level
#'   con <- c(0,0,1,-1,0,0)
#'
#' # checking estimability
#'   estimable(dmat, con)
#'
#'
estimable <- function(dmat, cmat) {
   res <- varFac(dmat, cmat, 1.e-5, "e")
   if (is.na(res)) return(FALSE)
   else return(TRUE)
}

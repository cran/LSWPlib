#' Wavelet Packet Basis for a single scale
#'
#' \code{get.flat.basis} returns the full set of packet indices relative to a basis
#' for a single scale from the wavelet packet table.
#'
#' @param scale The scale for which the indices of wavelet packet basis are returned.
#' Typically this is a positive integer.
#'
#' @details This function is used internally by other routines but it might be useful when
#' the wavelet packet spectral estimation over a fixed scale is of interest. The function returns
#' an object of class \code{lswpbb}, whose structure is the same to the object produced by \code{\link{best.basis}}.
#'
#' @return A matrix of two columns where each row refers to a different selected packet.
#' The first index is the argument \code{scale}, the second index refers to the packet within this level.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{best.basis}}, \code{\link{get.wavelet.basis}}
#'
#' @examples
#' get.flat.basis(scale = 4)
#'
#' @export

get.flat.basis <- function(scale)
{
  ans <- matrix(scale, nrow = 2^scale, ncol = 2)
  ans[,2] <- (1:(2^scale)) - 1
  class(ans) <- 'LSWPbasis'

  ans
}



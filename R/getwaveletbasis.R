#' Wavelet Packet indices for a Wavelet Basis
#'
#' \code{get.wavelet.basis} returns the full set of packet indices relative to a wavelet basis
#' selected from the wavelet packet table.
#'
#' @param scale The maximum scale for which the indices of the wavelet basis are returned.
#' Tipically this is a positive integer.
#'
#' @details This function is used internally by other routines but it might be useful when
#' the 'classical' wavelet spectral estimation is of interest. The function returns
#' an object of class \code{LSWPbasis}, whose structure is the same to the object produced by
#' \code{\link{best.basis}}.
#'
#' @return A matrix of two columns where each row refers to a different selected packet.
#' The first index refers to each given scale, the second index refers to the packet
#' within this level.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{best.basis}}, \code{\link{get.flat.basis}}
#'
#' @examples
#' get.wavelet.basis(scale = 4)
#'
#' @export

get.wavelet.basis <- function(scale)
{
  ans <- matrix(c(1:scale, scale), nrow = scale+1, ncol = 2)
  ans[,2] <- c(rep(2,scale),1) - 1
  class(ans) <- 'LSWPbasis'

  ans
}



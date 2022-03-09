#' Auto Convolution
#'
#'
#' \code{autoconv} computes the linear convolution of a numeric vector with itself.
#' It is based on the  \code{fft} function and is twicked to achieve maximum performance.
#'
#' @param x a real or complex vector
#'
#' @details The speed of calculation for the linear convolution depends upon the number of factors
#'          in the number representing the vector length. This implementation maximizes calculation speed for vectors
#'          of dyadic length, or lengths with a single factor.
#'
#' @return The linear auto convolution of a given vector with itself,
#'         which is equivalent with its inner product.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{fft}}, \code{\link{convolve}}.
#'
#' @examples
#'
#'  v <- rnorm(n = 64)
#' vv <- autoconv(x = v)
#'
#' @export

autoconv <- function(x)
{
  FFt <- stats::fft
  ny <- n1 <- length(x)
  y <- c(rep.int(0, n1), x)
  n <- length(x <- c(x, rep.int(0, n1)))
  z <- FFt(FFt(x) * Conj(FFt(y)), inverse = TRUE)
  ans <- Re(z)[-1]/n

  ans
}




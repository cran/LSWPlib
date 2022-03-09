#' Estimate an LSWP basis by penalised least squares
#'
#' \code{LSWPbasis} returns a matrix containing the wavelet packet basis indices.
#'
#' This function fits a wavelet packet basis to data using a penalised least square method.
#'
#' @param x a (locally stationary) time series of dyadic length.
#' @param wavelet the wavelet used to estimate the wavelet packet spectra.
#' @param lev.max the maximum scale for which the basis is fitted.
#' @param smooth should the penasised least squares cost functionals be smoothed? Default value is \code{TRUE}.
#' @param spa parameter for the local polynomial smothing implemented through \code{lowess}
#' @param penalty implemets increasing penalty for increasing scales.
#'
#' @details This function implements a data-driven basis selection of locally stationary time series.
#' The wavelet argument is specified as in other functions of this package.
#' Therefore, the current implementation allows for three discrete wavelets: Haar (\code{"haar"}),
#' Daubechies Extremal Phase linear filters of length 4 (\code{"d4"}) and Least Asymmetric linear filters of length 8 (\code{"la8"}).
#' Smoothing is controlled through the argument spa.
#'
#' @return A matrix of dimensions |b| x 2, where |b| is the number of packets in the basis.
#' The first column contains the scale indices of each packet in the basis, the second column contains the packet index within each scale.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{LSWPspec}}, \code{\link{LSWPsim}}
#'
#' @examples
#'
#' wpb <- LSWPbasis(x = sp500, wavelet = 'la8', lev.max = 4, smooth = TRUE, spa = 0.35)
#'
#' @export


LSWPbasis <- function(x, wavelet, lev.max, smooth, spa, penalty = 0.976)
{
  lx <- length(x)
  Jx <- floor(log(lx,2))
  wpct <- wpct.fin <- vector(mode = 'list', length = lev.max)
  stopifnot(lev.max <= lx)
  wpen <- rep(0, lev.max)

  for(lev in 1 : lev.max)
  {
    basis.lev <- get.flat.basis(scale = lev)
    AAb <- .Abb(bb = basis.lev, wavelet = wavelet, inverse = FALSE)
  #  AAAb <- abs(solve(AAb))
    ewper.lev <- LSWPspec(x, lev = lev, bb = basis.lev, wavelet = wavelet, smooth = smooth, spa = spa, correct = TRUE, AA = AAb)
    meanper.lev <- apply(ewper.lev, 1, mean, na.rm = TRUE)
    wpct[[lev]] <- -meanper.lev
    wpen[lev] <- -sum(meanper.lev)
  }

  wp.weights <- (sum(wpen) - wpen)/sum(wpen)

  for(lev in 1:lev.max) wpct.fin[[lev]] <- wpct[[lev]] * wp.weights[lev] * penalty^lev

  wp.basis <- best.basis(wpc = wpct.fin)
  class(wp.basis) <- 'LSWPbasis'

  wp.basis
}


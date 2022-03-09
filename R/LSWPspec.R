#' Locally Stationary Wavelet Packet Spectral Estimation
#'
#' \code{LSWPspec} returns the spectral estimate of a locally stationary time series
#' characterized by a wavelet packet basis.
#'
#' @param x a real valued numeric vector containing a time series of dyadic length.
#' @param lev the maximum level for which the spectra should be estimated.
#' @param bb a wavelet packet basis for which the spectra is estimated.
#' @param wavelet wavelet used to estimate the wavelet packet spectra. Possible values are \code{"haar"}, \code{"d4"} and \code{"la8"}. See also Details.
#' @param smooth logical. If \code{FALSE} the returned spectral estimate is not smoothed. Default value is \code{FALSE}. See also Details.
#' @param spa window length for spectral smoothing. Increasing values increase the smoothing.
#' @param correct logical. Should the returned spectral estimate be unbiased? Default is \code{TRUE}.
#' @param AA this argument is for internal use only and should be left alone. See also Details.
#'
#' @details The current implementation allow the use of these three well known Daubechies discrete wavelets for spectral estimation. Default choice is the
#'          \code{"la8"} wavelet which has decent control over frequency leakage characterizing compactly supported filters. In this initial implementation
#'          smoothing is provided by local polynomials through the \code{lowess} function and the smoothing parameter \code{spa} is passed to \code{lowess}.
#'          Future package versions will allow for different smoothing methods. The argument \code{AA} is tipically used by other functions to provide the
#'          inner product matrix when running simulations. For a direct usage on a single time series the matrix is calculated internally usig the default settings.
#'
#' @return A matrix containing the time-frequency spectral estimate where each column corresponds to a different time point and ech row
#' corresponds to a different packet from the given basis.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{get.wavelet.basis}}, \code{\link{LSWPsim}}, \code{\link{best.basis}}.
#'
#' @examples
#'
#'  wb <- get.wavelet.basis(4)
#' wpp <- LSWPspec(x = sp500, lev = 4, bb = wb, wavelet = 'la8', smooth = TRUE, spa = 0.35)
#'
#' @export

LSWPspec <- function(x, lev, bb, wavelet, smooth, spa, correct = TRUE, AA = NULL)
{
  if (wavelet == "haar") { fam = "DaubExPhase" ; num = 1 }
  if (wavelet == "d4") { fam = "DaubExPhase" ; num = 2 }
  if (wavelet == "la8") { fam = "DaubLeAsymm" ; num = 8 }

  if(ncol(bb) != 2) bb <- t(as.matrix(bb))

  nrb <- nrow(bb)
  lex <- length(x)
  indt <- 1:lex
  a1 <- wavethresh::wpst(x, filter.number = num, family = fam)
  a2 <- floor(log(lex, 2))
  a3 <- a5 <- matrix(nrow = nrb, ncol = lex)

  for(p in 1: nrb)
  {
    j <- bb[p,1]
    b <- bb[p,2]
    a3[p,] <- wavethresh::accessD.wpst(a1, level = ceiling(a2 - j), index = b)
  }
  a3 <- a3^2

  if(smooth == TRUE)
    for(i in 1 : nrb) a5[i,] <- stats::loess(formula = a3[i,] ~ as.numeric(indt), span = spa)$fitted
  else a5 = a3

  if(correct == TRUE)
  {
    if(is.null(AA)) AA <- .Abb(bb, wavelet = wavelet)
    Ai <- solve(AA)
    ans <- Ai %*% a5
  }

  else ans <- a5
  class(ans) <- 'LSWPspec'

  ans
}





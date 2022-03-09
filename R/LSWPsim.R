#' Simulation of LSWP processes
#'
#' \code{LSWPsim} returns simulated time series from a specified LSWP specification.
#'
#' This function produces one or multiple realizations of an LSWP process that is
#' specified in terms of a wavelet packet basis and its corresponding spectra.
#'
#' @param bb a wavelet packet basis for the simulated process.
#' @param spec a (locally stationary) spectra corresponding to the wavelet packet basis.
#' @param lev the maximum level of the basis that is considered for simulation. Usually this is set as the maximum level in \code{bb}.
#' @param wavelet the Daubechies wavelet used to build wavelet packets to simulate the process. See also Details.
#' @param N the number of realizations to be simulated.
#'
#' @details The function simulates realizations accordingly to the specified arguments.
#' The wavelet argument is specified as in other functions of this package.
#' Therefore, the current implementation allows for three discrete wavelets: Haar (\code{"haar"}),
#' Daubechies Extremal Phase linear filters of length 4 (\code{"d4"}) and Least Asymmetric linear filters of length 8 (\code{"la8"}).
#'
#' @return If \code{N = 1} the function returns a vector containing the simulated time series.
#' If \code{N > 1} the function returns a matrix with \code{N} columns each containing a different
#' simulated series.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{LSWPspec}}, \code{\link{best.basis}}, \code{\link{get.flat.basis}}, \code{\link{get.wavelet.basis}}.
#'
#' @examples
#'
#' wwb <- get.flat.basis(scale = 4)
#' wwp <- matrix(2^{-(1:4)}, nrow = 4, ncol = 512, byrow = FALSE)
#' xt <- LSWPsim(bb = wwb, spec = wwp, lev = 4, wavelet = 'la8', N = 10)
#'
#' @export

LSWPsim <- function(bb, spec, lev, wavelet, N)     # done
{

  if(any(spec < 0)) spec <- .resc.spec(spec)
  TT <- ncol(spec)
  rr <- nrow(spec)
  WP <- array( dim  = c(TT, rr, N))
  e <- specs <- array(stats::rnorm(TT*rr*N), dim = c(TT, rr, N))
  for(i in 1:N) specs[,,i] <- sqrt(t(spec)) * e[,,i]
  coe <- .wp.filt.seq(lev)

  for(i in 1:rr)
  {
    W <- matrix(0, ncol = TT, nrow = TT)
    fseq <- .pk.extract(r = i, bb = bb, J = lev, coe = coe)
    wfil <- .wavelet.filter2(wf.name = wavelet, filter.seq = fseq)
    lfi <- length(wfil);

    print(lfi)

    wfill <- c(wfil, rep(0, TT - lfi))

    for(k in 1:TT) W[k,] <- .circular(wfill, k-1)

    for(j in 1:N) WP[,i,j] <- W %*% specs[,i,j]
  }
         ans <- apply(WP, c(1,3), sum)
  class(ans) <- 'LSWPsim'

  ans
}





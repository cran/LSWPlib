#' Plot for Objects of Class LSWPspec
#'
#' \code{plot.LSWPspec} returns the plot for objects of class \code{"LSWPspec"}, typically (but not exclusively) a wavelet packet spectral estimate.
#'
#' @param x an object of class lswpspec.
#' @param y not used, is set to \code{NULL}.
#' @param ... not currently used.
#'
#' @details This function implements the \code{plot} method for objects of class \code{"LSWPspec"}. It is mainly used to plot spectral estimates returned by \code{LSWPspec}.
#' The label of the vertical axis uses the wavelet packet basis index notation \code{p = 1,2,...,|b|}, where |b| is the number of packets in a wavelet
#' packet basis as defined in Cardinali and Nason (2017). The label of the horizontal axis is the time index.
#'
#' @return A plot of the time-varying spectral estimates.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{LSWPspec}}, \code{\link{LSWPbasis}}
#'
#' @examples
#'
#'  wb <- get.wavelet.basis(4)
#' wpp <- LSWPspec(x = sp500, lev = 4, bb = wb, wavelet = 'la8', smooth = TRUE, spa = 0.35)
#'        plot(wpp)
#'
#' @export

plot.LSWPspec <- function(x, y, ...)
{
  mat <- as.matrix(t(x))
  bsize <- ncol(mat)
  colnames(mat) <- 1:bsize

  stats::plot.ts(mat, type = 'h', ann = TRUE, ylab = 'p', xlab = 't',
          nc = 1, main = '', sub = '', axes = FALSE, frame.plot = TRUE,
          mar.multi = c(0, 0, 0, 2.1),  oma.multi = c(6, 5, 5, 0), mgp = c(0,0,0))

}



#' Wavelet Packet Marginal Spectra
#'
#' \code{marg.spec} returns the time-average spectra of LSWP processes for each packet.
#'
#' This function computes the frequency intervals corresponding to each packet, along with the (time) average spectra for each packet.
#'
#' @param bas is a (|b| x 2) matrix containing indices of a wavelet packet basis, where |b| is the number of packets in the basis.
#' @param spec is a (|b| x T) matrix containing, in each row, the values of the time-varying spectra for each packet.
#' @param plot should a plot of the marginal spectra vs frequency intervals be returned?
#'
#' @details This function is used to compute, and eventually plot, the time averaged spectra (or spectral estimate) vs packet frequencies.
#' The arguments bas and spec shuld be provided as matrices.
#'
#' @return A (|b| x 2) matrix. In the first column the lower frequency relative to each packet is displayed. The second column contains the (time) average spectra.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{best.basis}}, \code{\link{get.wavelet.basis}}
#'
#' @examples
#'
#' wwb <- get.wavelet.basis(scale = 4)
#' wwp <- matrix(2^{-(1:4)}, nrow = 4, ncol = 512, byrow = FALSE)
#' msp <- marg.spec(bas = wwb, spec = wwp, plot = TRUE)
#'
#' @export

marg.spec <- function(bas, spec, plot = TRUE)
{
  fre <- function(bas) bas[2]/(2^(bas[1] + 1))

  if(is.matrix(spec))  msp <- apply(spec, 1, mean, na.rm = TRUE)
  if(!is.matrix(spec)) msp <- spec

  ome <- apply(bas, 1, fre)
  tmp <- cbind(ome, msp)
  ans <- as.data.frame(tmp[order(tmp[,1]),])
  nr <- nrow(bas)
  ome <- c(0, ans[,1])
  msp <- c(ans[,2], msp[nr])

  if(plot)
  {
    plot(ome, msp, type = 's', ylim = c(0, max(spec, na.rm = TRUE)), xlab = 'Frequency', ylab = 'Marginal Spectra')

    for(i in 1:length(ome)) graphics::abline(v = ome[i], lty = 2, col = 2)
  }

  anss <- cbind(ome, msp)

  return(anss)
}




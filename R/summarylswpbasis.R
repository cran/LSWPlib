#' Summary method for Objects of Class LSWPbasis
#'
#' \code{summary.LSWPbasis} returns a table containing the wavelet packet basis and its wavelet packet bases index notation.
#'
#' @param object an object of class \code{"LSWPbasis"}, typically (but not exclusively) returned by the function \code{LSWPbasis}.
#' @param ... not currently used.
#'
#' @details This function is used to print a wavelet packet basis with the wavelet packet basis index notation \code{p = 1,2,...,|b|},
#' where |b| is the number of packets in a wavelet packet basis as defined in Cardinali and Nason (2017). The doublets \code{"j_p, i_p"} refer,
#' to the scale and packet index within each scale, respectively.
#'
#' @return Print an object of class \code{LSWPbasis} with LSWP basis notation.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{LSWPbasis}}, \code{\link{best.basis}}
#'
#' @examples
#'
#' wpb <- LSWPbasis(x = sp500, wavelet = 'la8', lev.max = 4, smooth = TRUE, spa = 0.35)
#' summary(wpb)
#'
#' @export

summary.LSWPbasis <- function(object, ...)
{
  if(!(class(object) == 'LSWPbasis')) stop(paste(object, 'must be of class LSWPbasis'))

  nrb <- nrow(object)
  p <- 1:nrb
  ans <- data.frame(cbind(p,object))
  colnames(ans) <- c('p', 'j_p', 'i_p')

  cat(paste0('This is a LSWP `basis` of dimension |b| = ',nrb), "\n", fill = 1)

  print(ans, row.names = F)
}


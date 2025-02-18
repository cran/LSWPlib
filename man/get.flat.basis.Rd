% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getflatbasis.R
\name{get.flat.basis}
\alias{get.flat.basis}
\title{Wavelet Packet Basis for a single scale}
\usage{
get.flat.basis(scale)
}
\arguments{
\item{scale}{The scale for which the indices of wavelet packet basis are returned.
Typically this is a positive integer.}
}
\value{
A matrix of two columns where each row refers to a different selected packet.
The first index is the argument \code{scale}, the second index refers to the packet within this level.
}
\description{
\code{get.flat.basis} returns the full set of packet indices relative to a basis
for a single scale from the wavelet packet table.
}
\details{
This function is used internally by other routines but it might be useful when
the wavelet packet spectral estimation over a fixed scale is of interest. The function returns
an object of class \code{lswpbb}, whose structure is the same to the object produced by \code{\link{best.basis}}.
}
\examples{
get.flat.basis(scale = 4)

}
\references{
A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
            Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
}
\seealso{
\code{\link{best.basis}}, \code{\link{get.wavelet.basis}}
}
\author{
Alessandro Cardinali
}

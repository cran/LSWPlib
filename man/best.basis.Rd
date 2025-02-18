% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bestbasis.R
\name{best.basis}
\alias{best.basis}
\title{Best basis selection from a dyadic tree}
\usage{
best.basis(wpc)
}
\arguments{
\item{wpc}{this is a list containing the cost functionals associated with each packet. Each element in the list corresponds to a level of the dyadic tree.}
}
\value{
A matrix of two columns where each row refers to a different selected packet. The first index refers to the tree level, the second index refers to the packet within that level.
}
\description{
\code{best.basis} returns a selection of packets from a dyadic tree where the selection is made by minimizing the cost functionals associated with each packet.
}
\details{
The function implements best basis selection from an arbitrary tree. Typically this tree is produced by other functions in this package and this function
is also used to produce a selected basis. Since best basis methods are of general interest this function has been exported
for possible other uses.
}
\examples{

costs <- vector(mode = 'list', length = 4)
for(i in 1:4) costs[[i]] <- rnorm(2^i)^2
best.basis(wpc = costs)

}
\references{
A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
            Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
}
\seealso{
\code{\link{LSWPbasis}}, \code{\link{LSWPspec}}.
}
\author{
Alessandro Cardinali
}

#' Best basis selection from a dyadic tree
#'
#' \code{best.basis} returns a selection of packets from a dyadic tree where the selection is made by minimizing the cost functionals associated with each packet.
#'
#' @param wpc this is a list containing the cost functionals associated with each packet. Each element in the list corresponds to a level of the dyadic tree.
#'
#' @details The function implements best basis selection from an arbitrary tree. Typically this tree is produced by other functions in this package and this function
#' is also used to produce a selected basis. Since best basis methods are of general interest this function has been exported
#' for possible other uses.
#'
#' @return A matrix of two columns where each row refers to a different selected packet. The first index refers to the tree level, the second index refers to the packet within that level.
#'
#' @author Alessandro Cardinali
#'
#' @references A. Cardinali and G.P. Nason (2017). Locally Stationary Wavelet Packet Processes:
#'             Basis Selection and Model Fitting. Journal of Time Series Analysis, 38:2, 151-174.
#'
#' @seealso \code{\link{LSWPbasis}}, \code{\link{LSWPspec}}.
#'
#' @examples
#'
#' costs <- vector(mode = 'list', length = 4)
#' for(i in 1:4) costs[[i]] <- rnorm(2^i)^2
#' best.basis(wpc = costs)
#'
#' @export

best.basis <- function(wpc)
{
      wpc2 <- wpc
       lev <- length(wpc)
        bb <- matrix(nrow = 2^lev, ncol = 2)
        nb <- matrix(0, nrow = lev, ncol = 2^lev)
     dimme <- 1:2^lev
       ble <- 2^(1:lev)
         I <- .indix(lev)

for(j in (lev - 1) : 1){
              for(b in 1 : 2^j){

if(wpc[[j]][b] <= sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b]) && wpc2[[j]][b] <= sum(wpc2[[j+1]][(2*b) - 1], wpc2[[j+1]][2*b])) { nb[j, I[j,] == b] <- 1}


else if(wpc[[j]][b] > sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b]) && wpc2[[j]][b] > sum(wpc2[[j+1]][(2*b) - 1], wpc2[[j+1]][2*b]) )
                                         {
                                           wpc[[j]][b] <- sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b])

       if(sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b]) == sum(wpc2[[j+1]][(2*b) - 1], wpc2[[j+1]][2*b])) {nb[j+1, I[j+1,] == (2*b) | I[j+1,]  == (2*b) - 1] <- 1}
                                         }

else if(wpc[[j]][b] <= sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b]) && wpc2[[j]][b] > sum(wpc2[[j+1]][(2*b) - 1], wpc2[[j+1]][2*b]))
                                         {
                                           wpc[[j]][b] <- sum(wpc[[j+1]][(2*b) - 1], wpc[[j+1]][2*b])
                                           nb[j+1, I[j+1,] == (2*b) | I[j+1,]  == (2*b) - 1] <- 0
                                         }
                               }
                       }

for(b in 1 : (2^lev)){
           for(j in (1 : (lev - 1))){
                                      if(nb[j,b] == 1){ nb[min(lev,(j+1)):lev, b] <- 0; j <- 1; b <- min(b+1, 2^lev)}
                                    }
                     }

pl <- pc <- rep(0, lev+1)
pc[1] <- pl[1] <- 1

for(j in 2 : (lev+1)){
jj <- j - 1

                   po <- seq(1, 2^lev, 2^(lev-jj))
                   pa <- I[jj, nb[jj,] == 1]

 if(length(pa) == 0) j <- j+1

else                              {
                      pu <- as.vector(stats::na.omit(pa[po]))
                   pl[j] <- length(pu)
                      pc <- cumsum(pl)
   bb[pc[jj]:(pc[j]-1),] <- cbind(rep(jj, pl[j]), pu)


                                  }
                     }

          bb <- tmp <- matrix(stats::na.omit(bb), ncol = 2)
      bb[,2] <- tmp[,2] - 1
dimnames(bb)[[2]] <- c("scale", "packet")

bb
}






# internals for LSWPlib package




.convert.bb <- function(bb)
{
       ans <- bb
   ans[,2] <- bb[,2] + 1

  ans
}




#




.wp.filt.seq <- function(J)
{
# sequenze of highpass (1) and lowpass (0) filters for a wavelet packet table #
# the order of rows is inverse because so is required by "wavelet.filters2"#

        n <- 2^J
    matri <- matrix(rep(0, (J * (2^J))), nrow = J)

for(j in 1 : J)
      {
              m <- 2^j
             mn <- n / m
              s <- seq(1:m)
              b <- seq(1, n, mn)

for(k in 1:m) {if(.mod(s[k],4) == 2 | .mod(s[k],4) == 3) matri[j,c(b[k]:(b[k]+(mn -1)))] <- rep(1,mn)}
      }

matri[J:1,]
}





#




.pk.extract <- function(r, bb, J, coe)
{
coe <- as.matrix(coe)

if(ncol(bb) != 2) bb <- t(as.matrix(bb))
if(ncol(coe) < 2) coe <- t(coe)

bb <- .convert.bb(bb)
a1 <- coe[J:1,]
a2 <- bb[r,]
 j <- a2[1]
 b <- a2[2]
ans <- (coe[1:j, b])

ans
}




#




.wavelet.filter2 <- function(wf.name, filter.seq)
{
    cascade <- function(f, x, j)
    {
        L <- length(f)
        N <- length(x)
        M <- (L - 1) * 2^j
        M1 <- M - L + 2
        M2 <- 2 * M - L + 2
        if (N > M1)
            stop("x is too long\n")
        else x <- c(x, rep(0, M1 - N))
        xj <- c(rep(0, M), x, rep(0, M))
        yj <- rep(0, M2)
        for (i in 1:L) yj <- yj + f[L - i + 1] * xj[1:M2 + (i - 1) * 2^j]
        yj
    }
    if (is.character(wf.name))
        wf <- waveslim::wave.filter(wf.name)
    else wf <- wf.name
    J <- length(filter.seq)
    key <- rev(filter.seq)     # was rev() #
    f <- 1
    fl <- wf$lpf
    fh <- wf$hpf
    for (k in 1:J)
    {
          if (key[k] == 1) f <- cascade(fh, f, k - 1)
        else if (key[k] == 0) f <- cascade(fl, f, k - 1)
        else stop("Invalid filter.seq")
    }

    f
}






#





.acwp <- function(r1, bb, wavelet)
{
      mbb <- J <- max(bb[,1])
      coe <- .wp.filt.seq(J)
     lwav <- switch(wavelet, 'haar' = 2, 'd4' = 4, 'la8' = 8)
     fse1 <- .pk.extract(r1, bb, J, coe)
      hb1 <- .wavelet.filter2(wavelet, filter.seq = fse1)
      lb1 <- length(hb1)
    ac.wp <- autoconv(x = hb1)
     lfil <- floor((length(ac.wp) - 1)/2)
     mfil <- ((2^mbb - 1)*(lwav - 1) + 1)*2 - 1
    centr <- (mfil + 1)/2
      ans <- rep(0, mfil)

            ans[(centr-lfil):(centr+lfil)] <- ac.wp

  ans
}






#





.wpfilt <- function(r1, bb, wavelet, stepf = F, steps)
{
      mbb <- J <- max(bb[,1])
      coe <- .wp.filt.seq(J)
     lwav <- switch(wavelet, 'haar' = 2, 'd4' = 4, 'la8' = 8)
     fse1 <- .pk.extract(r1, bb, J, coe)
      hb1 <- .wavelet.filter2(wavelet, filter.seq = fse1)
      ans <- c(0,hb1,0)

      if(stepf == T) ans <- stats::stepfun(x = steps, y = c(0,hb1,0))

ans
}




#----------






.Aacwp = function(bb, wavelet)
{
           rr <- nrow(bb)
         rind <- 1:rr
     acwp.mat <- sapply(X = rind, FUN = .acwp, bb = bb, wavelet = wavelet)
          ans <- t(acwp.mat) %*% acwp.mat

ans
}




#----------




.Abb <- function(bb, wavelet, inverse = FALSE)
{
        ans <- A <- .Aacwp(bb = bb, wavelet = wavelet)

        if(inverse == T) ans <- solve(A)

ans
}




#---------




.mod <- function(a,b)
{
q <- floor(a/b)
r <- a - q*b

r
}




#---------




.circular <- function(v, i)
{
if(i == 0) w <- v
else
{
            n <- length(v)
            w <- rep(0, n)
   w[(1+i):n] <- v[1:(n-i)]
       w[1:i] <- v[(n - i + 1):n]
}

w
}




#




.resc.spec <- function(ewps)
{
   n.lev <- nrow(ewps)
      TT <- ncol(ewps)
      tt <- 1:TT
       a <- b <- rep(0, n.lev)
       a <- apply(ewps, 1, stats::var)
           ews[ewps < 0 | is.na(ewps)] <- 0
       b <- apply(ewps, 1, stats::var)
      re <- a/b
           if(any(is.na(re))) re[is.na(re)] <- 1
           if(any(re == Inf)) re[re == Inf] <- 1
    ans <- ewps * sqrt(re)

ans
}




#




.indix <- function(J)
{
ans <- matrix(nrow = J, ncol = 2^J)
a2 <- 1:J

for(j in 1:J)
         {
             a3 <- 1 : 2^j
             a4 <- length(a3)
             ans[j,] <- rep(a3, rep(2^(J - j), a4))
         }
ans
}





#





.wp.fit <- function(bas, wpt)
{
wp.select <- function(x,y){y[[x[1]]][x[2]]}
      bas <- .convert.bb(bas)
      wpb <- apply(bas, 1, wp.select, wpt)
      ans <- sum(wpb, na.rm = TRUE)

c(ans, wpb)
}

































































































































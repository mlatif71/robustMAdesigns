#' @title Estimate and test a contrast matrix/vector
#'
#' @description Estimate and test a contrast matrix/vector
#'
#' @param z vector of responses
#' @param dmat design matrix of interest
#' @param cmat contrast matrix corresponding to the effect of interest
#' @param alternative specification of alternative hypothesis
#'
#' @return Returns a vector of test statistics, estimate of the contrast  matrix/vector, p--value, mean squared error, etc. Each row of the design  matrix corresponds to a response
#'
#' @export contrastEst
#'
#' @examples
#' # design matrix
#'   dmat <- desMatrix("1x3", design="CL")
#'   dmat <- rbind(dmat, dmat)
#'
#' # contraxt vector for comparing first and second factor level
#'   con <- c(0,0,1,-1,0)
#'
#' # response
#'   y <- rnorm(6)
#'   contrastEst(z=y, dmat=dmat, cmat=con)
#'
#' @seealso \code{\link{desMatrix}}, \code{\link{contMatrix}}
#'
#' @references Searle SR (1971). Linear Models. Wiley.
#'
contrastEst <-  function(z, dmat, cmat, alternative="two.sided") {
     x <- dmat
     cc <- cmat
     z <- matrix(z, nrow = length(z), ncol = 1)
     if (is.null(dim(cc))) {
       len <- 1
     } else {
       len <- dim(cc)[1]
     }
     x <- as.matrix(x)
     xt <- t(x)
     n <- nrow(x)
     xx <- xt %*% x
     ff <- n - qr(x)$rank
     if (ff < 1) {
       message("Error df is too small!")
       return(NA)
     }
     gxx <- MASS::ginv(xx, tol = sqrt(.Machine$double.eps)/10)
     sse <- t(z) %*% (diag(n) - x %*% gxx %*% xt) %*% z
     mse <- sse/ff
     if (len > 1) {
        if (!estimable(x, cc)) {
          message("Effect is not estimable")
          return(NA)
        }
        cc <-  as.matrix(cc)
        cb <- cc %*% gxx %*% xt %*% z
        cxxc <- cc %*% gxx %*% t(cc)
        ff1 <- qr(cxxc)$rank
        test.stat <- t(cb) %*% MASS::ginv(
          cxxc,
          tol = sqrt(.Machine$double.eps)/10) %*% (cb) / (mse * ff1)
        pval <- 1 - stats::pf(test.stat, ff1, ff)
        df1 <- ff1
        df2 <- ff
        cb <- NA
        res <- round(c(test.stat, mse, df1, df2, pval),3)
        names(res) <- c("F", "mse", "df.num", "df.denom", "p.val")
        return(res)
     } else {
        if (!estimable(x, cc)) {
          message("Effect is not estimable!!")
          return(NA)
        }
        cb <- cc %*% gxx %*% xt %*% z
        cc <- as.vector(cc)
        test.stat <- cb/sqrt(mse*(cc %*% gxx %*% cc))
        if (alternative == "two.sided") {
          pval <- min(2*(1 - stats::pt(abs(test.stat), ff)), 1)
        } else if (alternative == "less") {
          pval <-  stats::pt(test.stat, ff)
        } else if (alternative == "greater") {
          pval <- 1 - stats::pt(test.stat, ff)
        }
        df1 <- ff
        df2 <- NA
        res <- c(round(cb,3), round(test.stat,3), round(mse,3), round(df1,3),
                 round(pval,3))
        names(res) <- c("cont.est", "t", "mse", "df", "p.val")
        return(res)
     }
     res <- c(cb, test.stat, mse, df1, df2, pval)
     return(res)
 }


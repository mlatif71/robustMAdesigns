comb <- function(n, r, v = 1:n) {
  #
    if (r == 1) {
      matrix(v, n, 1)
    } else if (r == n) {
      matrix(v, 1, n)
    } else {
      rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])),
               Recall(n - 1, r, v[-1]))
    }
}
#
varFac <- function(x, cc, tol=1.e-5, type=c("a", "e", "d", "t")) {
  #
  x <- as.matrix(x)
  # make sure contrast is a matrix
  if (is.null(dim(cc))) cc <-  t(as.matrix(cc))
  else cc <- as.matrix(cc)
  #
  if (is.null(dim(x))) xx <- outer(x, x)
  else xx <- t(x) %*% x # this is for designs with one slide
  gxx <- MASS::ginv(xx, tol = sqrt(.Machine$double.eps) / 10)
  # estimability
  tem <- cc %*% gxx %*% xx - cc
  # variance factor
  # dispersion matrix
  cxxc <- cc %*% gxx %*% t(cc)
  cxxc1 <- cc %*% xx %*% t(cc)

  # added qr() to compute rank of a matrix on 15.04.2004
   ri <- qr(cxxc)$rank
  #
  atol <- sum(diag(tem %*% t(tem)))
  #
  # norm cc' : c is row vector or matrix
  cc.norm <- sum(diag(cc %*% t(cc)))

  # eigen value of cxxc
  # without using symmetric=T sometime eigen values are comples for 3x2 design with 4 slides are missing
  # only.values=T only make the codes faster
  ei.cxxc <- eigen(cxxc, symmetric = TRUE,
                   only.values = TRUE)$values

  if (atol <= tol) {
     if (type == "a") return(ri / sum(diag(cxxc)))
     # minimax for dispersion matrix or maximin of information matrix
     else if (type == "e") return(cc.norm / max(ei.cxxc))
     else if (type == "d") {
         if (sum(ei.cxxc > 1.e-4)) {
           return(1 / prod(ei.cxxc[ei.cxxc > 1.e-4])^(1 / ri))
         } else {
           stop("None of the eigen value is positive")
         }
     } else if (type == "t") return(mean(diag(xx)))
  }
  else return(NA)
}

# centering matrix
 Pmat <- function(p) {
   diag(p) - (1 / p) * matrix(1, nrow = p, ncol = p)
   }
# sum vector
 one <- function(p) {
   matrix(1, nrow = p, ncol = 1)
   }
#
circleLoop <- function(n){
  out <- matrix(rep(1:n, each = 2)[-c(1, 2 * n)],
                ncol = 2, byrow = TRUE)
  out <- rbind(out, c(n, 1))
  return(out)
}

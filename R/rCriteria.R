#' @title robustness criteria
#'
#' @description Computes robustness criteria for a given design and contrast matrix/vector
#'
#"
#' @param dmat design matrix of interest
#' @param cmat contrast matrix corresponding to the effect of interest
#' @param cinfo a vector represents the number of rows corresponding to the
#' contrast matrices of interest
#' @param type specifies the type of the efficiency criteria, either
#' "a" or "e" or "d"
#' @param cname name of the contrasts
#'
#' @export rCriteria
#'
#' @return A list of \code{bdn}, \code{avgEff}, and \code{pED}

#'
#' @references
#' Latif AHMM, Bretz F, and Brunner E (2009). Robustness consideration
#'   in selecting efficient two-color microarray designs. Bioinformatics,
#'    25(18), 2355-2361.
#'
rCriteria <- function(dmat, cmat, cinfo=NULL, type="e",
                      cname=NULL) {
  #
   if (is.null(cinfo)) {
      cinfo <- ifelse(is.null(dim(cmat)), 1, dim(cmat)[1])
   }
  #
   nc <- length(cinfo)
   bdn <- rep(0, nc)
   #
   frun <-  unlist(eCriteria(dmat = dmat,
                             cmat = cmat,
                             m = 0,
                             cinfo = cinfo,
                             type = type)$eff)
   #
   avgEff <-  rbind(NULL, frun)
   pED <- rbind(NULL, !is.na(frun))
   ex.sl <- 1
   #
   while (sum(bdn == 0) > 0) {
      res <- eCriteria(dmat = dmat,
                       cmat = cmat,
                       m = ex.sl,
                       cinfo = cinfo,
                       type = type)$eff
      #
      nr.res <- ifelse(nc > 1,  nrow(res), length(res))
      #
      if (nc > 1) {
         avgEff <- rbind(avgEff,
                         apply(res, 2, sum, na.rm = TRUE) / nr.res)
         pED  <- rbind(pED,
                       apply(res, 2,
                             function(x) sum(!is.na(x))) / nr.res)
         res1 <- is.na(apply(res, 2, sum))
      } else {
         avgEff <- c(avgEff, sum(unlist(res), na.rm = T) / nr.res)
         pED <- c(pED, sum(!is.na(res)) / nr.res)
         res1 <- is.na(sum(res))
      }
      #
      bdn1 <- sapply(res1, function(x) ifelse(x == 1, ex.sl, 0))
      bdn[bdn == 0] <- bdn1[bdn == 0]
      ex.sl <- ex.sl + 1
   }
   bdn <- as.data.frame(t(bdn))
   row.names(bdn) <- " "
#   return(avgEff)
   avgEff <- as.data.frame(round(avgEff, 3))
   pED <- as.data.frame(round(pED, 3))
   row.names(avgEff) <- row.names(pED) <- 0:max(bdn)
   #
   if (is.null(cname)) {
     colnames(avgEff) <- colnames(pED) <-
     colnames(bdn) <- paste("C", 1:nc, sep = "")
   } else {
     colnames(avgEff) = colnames(pED) = colnames(bdn) = cname
   }
   return(list(bdn = bdn, avgE = avgEff, pED = pED))
}
#######

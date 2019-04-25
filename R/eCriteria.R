#' @title Efficiency Criteria
#'
#' @description Computes efficiency criteria for a given design and the contrast matrix/vector
#'
#'
#' @param type types of efficiency criteria, e.g. "e" is for E-efficiency, which is the only available option in the current version. Other options such as "d" (D-efficiency) will be added later.
#'
#' @param cname names of the contrasts
#'
#' @param  dmat design matrix of interest
#'
#' @param  cmat contrast matrix corresponding to the effect of intererst
#'
#' @param cinfo a vector, represents the number of rows corresponding to the
#'  contrast matrices of interest
#'
#' @param m number of missing values to be considered for computing
#'  efficiency  criteria
#'
#' @seealso \code{contMatrix}, \code{desMatrix}, \code{estimable}
#'
#' @return A \code{list} of efficiency criteria (\code{eff}) and corresponding
#' set of missing arrays (\code{m})
#'
#' @export eCriteria
#'
#' @examples
#' # design matrix corresponding to a circular loop design from 3x2 layout
#'   x <- desMatrix(layout = "3x2", design="CL")
#'
#' # contrast matrices for the main effects of two factors
#'    ccA <- contMatrix(layout="3x2", effect="mainA")
#'    ccB <- contMatrix(layout="3x2", effect="mainB")
#'
#' # combining the contarst matrices
#'    cc = rbind(ccA, ccB)
#'
#' # corresponding value of cinfo
#'   cin = c(nrow(ccA), nrow(ccB))
#'
#' # e-optimality without any missing value
#'   eCriteria(dmat=x, cmat=cc, cinfo=cin, m=0)
#'
#' # e-efficiency with at most two missing values
#'   eCriteria(dmat=x, cmat=cc, cinfo=cin, m=2)
#'
#' @references
#'   Pukelsheim F (1993). Optimal designs of experiments. Wiley.
#'
#'   Landgrebe J, Bretz F, Brunner E (2005). Efficient design and
#'   analysis of two-color factorial microarray design.
#'   Computational Statistics & Data Analysis, 50, 499-517.
#'
#'   Latif AHMM, Bretz F, and Brunner E (2009). Robustness consideration
#'   in selecting efficient two-color microarray designs. Bioinformatics,
#'    25(18), 2355-2361.
#'
eCriteria <- function(dmat, cmat, cinfo=NULL, m=0, type="e", cname=NULL) {
    tol <- 1.e-5
    design.list <- as.matrix(dmat)
    xsl <- m
    if (!is.null(dim(cmat))) {
       cmat <- as.matrix(cmat)
       if (is.null(cinfo)) cinfo <- dim(cmat)[1]
    }
    else {
       cmat <- as.vector(cmat)
       if (is.null(cinfo)) cinfo <- 1
    }
    # main loop
    if (xsl == 0) {
        res <- Fun1.eCriteria(
          design.list = design.list, cmat = cmat,
          cinfo = cinfo, tol = tol, type = type
          )
#        res1 = c(res)
        ncon <- length(cinfo)
        res1 <- as.data.frame(t(c(res)))
        if (is.null(cname)) colnames(res1) <- paste("C", 1:ncon, sep = "")
        else  colnames(res1) <- cname
        row.names(res1) <- " "
        return(list(eff = res1))
        #return(c(res1, xSlides=0))
    }
    else {
        nc <- length(cinfo)
        sl1 <- comb(nrow(design.list), xsl)
        res <- as.matrix(t(apply(sl1, 1, Fun1.eCriteria, design.list, cmat,
                                 cinfo, tol, type)))
        eff <- res[,1:nc]
        arr <- res[, (nc + 1):(nc + xsl)]
        return(list(eff = eff, marrays = arr))
    }
     return()
}

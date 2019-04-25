#' @title Genetic Algorithm for Microarray Designs
#'
#' @description Genetic Algorithms for selecting near-optimal microarray designs for a specific number of arrays. Instead of specifying experimetnal layout, the design matrix corresponding to  the list  of arrays from which near-optimal designs will be searched can also be used.
#'
#'
#' @references Latif AHM and Brunner E (2016). A genetic algorithm for designing microarray experiments. Computational Statistics, 31(2), 499-424.
#'
#' @param dmat design matrix
#' @param layout design layout
#' @param cmat contrast matrix/vector
#' @param cinfo dimension of contrast matrix/vector
#' @param n number of available arrays
#' @param Pcross crossover probablity
#' @param crossType crossover type, currently available options are: \code{"one.point"}, \code{"two.points"}, and \code{"uniform"}
#' @param  Pmut mutattion probablity
#' @param selectType selection type, currently available options are: sampling proportional to fitness (\code{"SPF"}) and reminder stochastic sampling (\code{"RSS"})
#' @param  verbose controls the print out during iterations
#' @param  wt a vector of weights to indicate individual importance of  effects of interest
#' @param popSize population size of a generation
#' @param maxGen maximum number of generations algorithm will run unless the convergence criteria is met
#' @param burnOut number of generations for burning out
#' @param convergeNum convergence criteria
#' @param nIter number of iterations
#' @param m number of missing values to be considered for computing  fitness fuction
#' @param keep.elite preserving the best design of each generation
#' @param scaling linear scaling on fitness
#' @param cname contrast names
#'
#' @importFrom stats mad median pf pt weighted.mean
#'
#' @export madesign_ga
#'
#' @examples
#' # near-optimal designs with 7 arrays for the layout 3x2 and main effects are of interest
#'
#' # contrast matrix
#'   A <- contMatrix(layout = "3x2", effect = "mainA")
#'   B <- contMatrix(layout = "3x2", effect = "mainB")
#'   cmat <- rbind(A, B)
#'   cinfo <- c(nrow(A), nrow(B))
#'
#' \dontrun{
#' madesign_ga(n=7, cmat=cmat, cinfo=cinfo, layout="3x2")
#' }
#'
  madesign_ga <- function(dmat=NULL, layout=NULL, n, cmat, cinfo=NULL, Pcross=.75,
            crossType="one.point", Pmut=NULL, selectType="RSS",
               verbose=FALSE, wt=NULL, popSize=50,  maxGen=2000, burnOut=100,
                   convergeNum=10, nIter=1, m=0, keep.elite=FALSE,
                        scaling=TRUE, cname=NULL){
#######
     narrays <- n
     nmarrays <- m

     ## defining weight function for averaging opt over questions
     if (is.null(wt))  leaveOut <- rep(1, length(cinfo))
     else leaveOut <- wt
#
     if (verbose & nIter > 1) {
       cat("\n     Iterations:", 1, " of ", nIter, "\n")
     }
     # keep best designs from each generation
     des <- list()
     # if mutation prob null
     if (is.null(Pmut)) Pmut <- 1 / narrays
  #
     out <- GA.default(dmat, layout, narrays, cmat, cinfo, leaveOut, popSize,
              Pcross, Pmut, maxGen, verbose, crossType, selectType, burnOut,
                 convergeNum, nmarrays, keep.elite, scaling, cname)
     des[[1]] <- out$best$design
     nGen <- out$history$nGen
     if (!is.null(cname)) names(out$best$opt) <- cname
  #
     if (nIter > 1) {
       	for (i in 2:nIter) {
       	  if (verbose) {
       	    cat("   Iterations:", i, " of ", nIter, "\n")
       	  }
       	  #
       	  out1 <- GA.default(dmat, layout, narrays, cmat, cinfo, leaveOut,
                       popSize, Pcross, Pmut, maxGen, verbose, crossType,
                       selectType, burnOut, convergeNum, nmarrays, keep.elite,
                               scaling, cname)
       	  des[[i]] <- out1$best$design
	        nGen <- c(nGen, out1$history$nGen)
	        if (out1$best$opt.overall > out$best$opt.overall) {
	          out <- out1
	          }
	        }
	        opt <- t(sapply(des,
	                        function(x) unlist(
	                          eCriteria(dmat = out$input$dmat[x,],
	                                    cmat = out$input$cmat,
	                                    cinfo = out$input$cinfo)$eff)))
	        opt.overall <- apply(opt, 1, stats::weighted.mean, leaveOut, na.rm = T)
          #opt = data.frame(opt)
	        if (!is.null(cname)) colnames(opt) <- cname
	        best <- list(design = do.call("rbind", des),
	                     opt = opt,
	                     opt.overall = opt.overall,
	                     nGen = nGen)
	       out$best <- best
	     }
       invisible(out)
  }

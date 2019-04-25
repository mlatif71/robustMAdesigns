#' @title Generating a contrast matrix or vector
#'
#' @description Generates a contrast matrix/vector corresponding to a specific effect
#'
#'
#'
#' @param layout experimetnal layout of interest, e.g. "1x3" is for one-factor experiment with three treatments and "2x2" is for two-factor experiments with both the factors with two levels
#'
#' @param effect  depending on the \code{layout}, \code{effect} could be different, e.g. for one-factor experiment, "all-pair" and "global" can be defined as \code{effect} and for two-factor experiment, effcet could be main effects ("A" or "B"), interaction ("AxB"), and simple effects ("simA" or "simB")

#' @param commRef indicates whether the associated design is common reference
#'
#' @return Contrast matrix/vector is returned, whose first two columns/elements correspond to the dye effects
#'
#' @examples
#'
#' # For a 1x4 layout, contrast matrix corresponds to all pairwise comparisons
#'    contMatrix(layout="1x4", effect="all-pair")
#'
#' # For a 4x3 layout, contrast matrix corresponding to interaction
#'    contMatrix(layout="4x3", effect="AxB")
#'
#' @export contMatrix
#'
#'
contMatrix <- function(layout, effect, commRef=FALSE) {
  type <- effect
  dname <- unlist(strsplit(layout, split = "x"))
  p <- as.numeric(dname[2])
  if (p < 2) {
    message("Layout is not properly specified!")
    return(NA)
  }
  if (dname[1] == 1) {
     if (type == "all-pair") {
        cmat <- desMatrix(layout, design = "DS")
        cmat <- cmat[1:(nrow(cmat)/2), ]
        if (commRef) {
          cmat <- data.frame(cbind(cmat,0))
          colnames(cmat)[3:(p + 3)] <- c(as.character(1:p),"R")
          }
        cmat[,1:2] <- 0
        return(cmat)
     } else if (type == "global") {
#        p = as.numeric(dname[2])
        cmat <- Pmat(p)
        if (commRef) {
           cmat <- as.data.frame(cbind(0,0,cmat*p,0))
           colnames(cmat)[-c(1,2)] <- c(as.character(1:p), "R")
         } else {
           cmat <- as.data.frame(cbind(0,0,cmat*p))
           colnames(cmat)[-c(1,2)] <- as.character(1:p)
        }
        colnames(cmat)[1:2] <- c("cy3", "cy5")
        return(cmat)
     } else {
       message("Contrast matrix for this effect is not available!")
       return(NA)
     }
  } else if (dname[1] > 1) {
     p1 <- as.numeric(dname[1])
     p2 <- as.numeric(dname[2])
     p12 <- p1*p2
     nam <- lapply(1:p1, function(i) paste(i, 1:p2, sep = ""))
     if (type == "AxB") {
       cmat <- (kronecker(Pmat(p1), Pmat(p2))*p12)[seq(1,p12, by = 2),]
     } else if (type == "mainA") {
       cmat <- kronecker(Pmat(p1), (1/p2)*t(one(p2)))*p12
     } else if (type == "mainB") {
       cmat <- kronecker((1/p1)*t(one(p1)), Pmat(p2))*p12
     } else if (type == "simA") {
       cmat0 <- kronecker(Pmat(p1), diag(p2))*p1
       cmat <- list()
       for (i in 1:p2) {
         if (commRef) {
           cmat[[i]] <- as.data.frame(cbind(0, 0,
             cmat0[seq(i,p12,by = p2),], 0))
           colnames(cmat[[i]])[-c(1,2)] <- c(unlist(nam),"R")
           } else {
             cmat[[i]] <- as.data.frame(cbind(0, 0,
               cmat0[seq(i,p12, by = p2),]))
             colnames(cmat[[i]])[-c(1,2)] <- unlist(nam)
           }
         colnames(cmat[[i]])[1:2] <- c("cy3", "cy5")
       }
       return(cmat)
     } else if (type == "simB") {
       cmat0 <- kronecker(diag(p1), Pmat(p2))*p2
       cmat <- list()
       for (i in 1:p1) {
         if (commRef) {
           cmat[[i]] <- as.data.frame(
             cbind(0,
               0,
               cmat0[(1 + (i - 1)*p2):(i*p2),],0))
           colnames(cmat[[i]])[ -c(1,2)] <- c(unlist(nam),"R")
         } else {
           cmat[[i]] <- as.data.frame(
             cbind(0,0,cmat0[(1 + (i - 1)*p2):(i*p2),]))
           colnames(cmat[[i]])[-c(1,2)] <- unlist(nam)
         }
         colnames(cmat[[i]])[1:2] <- c("cy3", "cy5")
       }
       return(cmat)
     } else {
       message("Contrast matrix for this effect is not available!")
       return(NA)
     }
     if (commRef) {
       cmat <- as.data.frame(cbind(0,0,cmat,0))
       colnames(cmat)[-c(1,2)] <- c(unlist(nam),"R")
     } else {
       cmat <- as.data.frame(cbind(0,0,cmat))
       colnames(cmat)[-c(1,2)] <- unlist(nam)
     }
     colnames(cmat)[1:2] <- c("cy3", "cy5")
     return(cmat)
  }
}



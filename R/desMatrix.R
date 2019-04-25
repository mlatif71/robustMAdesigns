#' @title Generating a Design Matrix
#'
#' @description Generates a design matrix corresponding to an experimental layout
#' and a design
#'
#'
#' @param layout experimetnal layout of interest, e.g. "1x3" is for one-factor experiment with three treatments and "2x2" is for two-factor experiments with both the factors with two levels
#'
#'
#' @param design design from the layout of interest, e.g. dye-swarp ("DS"), circular loop ("CL"), etc.
#'
#' @export desMatrix
#'
#' @return A data frame of which first two columnscorrespond to the dye effects cy3 and cy5, respectively and the remaining columns correspond to the number of treatments to be compared.
#'
#' @examples
#'
#' # design matrix for dye-swap design from 1x3 layout
#'     desMatrix(layout="1x3", design="DS")
#'
#' # design matrix for circular loop design from 2x2 layout
#'     desMatrix(layout="2x2", design="CL")
#'
desMatrix <- function(layout,  design="DS") {
  dname <- unlist(strsplit(layout, split = "x"))
  if (length(dname) < 2) {
    message("Input for layout ",layout, " is not of correct format!")
    return(NA)
  }
  # 1xp layout
  if (dname[1] == 1) {
    p <- as.numeric(dname[2])
    if (design == "CL") {
       ind <- circleLoop(p)
       arrays <- matrix(0, nrow = p, ncol = p)
       arrays[col(arrays) == ind[,1]] <- 1
       arrays[col(arrays) == ind[, 2]] <- -1
    }
    else if (design == "CR") {
       arrays <- as.data.frame(cbind(1, -1, diag(p), -1))
       colnames(arrays) <- c("cy3", "cy5", as.character(1:p), "R")
       return(arrays)
    }
    else if (design == "DS") {
      arrays <- matrix(0, nrow = choose(p,2), ncol = p)
      ind <- comb(n = p, r = 2)
      arrays[col(arrays) == ind[,1]] <- 1
      arrays[col(arrays) == ind[,2]] <- -1
    #if(dyeSwap)
      arrays <- rbind(arrays, -1 * arrays)
    }
    arrays <- as.data.frame(cbind(1, -1, arrays))
    colnames(arrays) <- c("cy3", "cy5", as.character(1:p))
    return(arrays)
  }
  # n x m layouts
  else if (dname[1] > 1) {
    p <- as.numeric(dname[1]) * as.numeric(dname[2])
    arrays <- desMatrix(paste(1, p, sep = "x"), design = "DS")
    nam <- lapply(1:as.numeric(dname[1]),
                  function(i) paste(i, 1:as.numeric(dname[2]), sep = ""))
    #
    if (design == "CR") {
       arrays <- as.data.frame(cbind(1,-1,diag(p),-1))
       colnames(arrays) <- c("cy3", "cy5", unlist(nam), "R")
       return(arrays)
    }
    else if (design == "DS") {
       arrays <- arrays
    }
    else if (layout == "2x2") {
       if (design == "CL") arrays <- arrays[c(1,5,12,8),]
       else if (design == "XS") arrays <- arrays[c(3,9,4,10),]
       else if (design == "AS") arrays <- arrays[c(2, 8, 5, 11),]
       else if (design == "BS") arrays <- arrays[c(1,7,6,12),]
       else if (design == "XL") arrays <- arrays[c(1,4,6,9),]
       else {
         message("Design matrix for the design ",
                 design, " yet to be implemented!")
         return(NA)
       }
    }
    else if (layout == "3x2") {
       if (design == "CL") arrays <- arrays[c(1,7,14,30,26,17),]
       else if (design == "BS") arrays <- arrays[c(1,16,10,25,15,30),]
       else if (design == "AL") arrays <- arrays[ c(4, 26, 17, 9, 29, 22),]
       else if (design == "XL") arrays <- arrays[ c(3, 13,23, 6,12,20),]
       else if (design == "TL") arrays <- arrays[c(4,28,18,6,12,24),]
       else if (design == "RS") arrays <- arrays[c(5,20, 8,23,10,25),]
       else {
         message("Design matrix for the design ",
                 design, " yet to be implemented!")
         return(NA)
       }
    }
    else {
      message("Design matrix for ", design, " of ",
              layout, " yet to be coded!")
      return(NA)
    }
    #
    arrays <- as.data.frame(arrays)
    colnames(arrays)[-c(1,2)] <- unlist(nam)
    colnames(arrays)[1:2] <- c("cy3", "cy5")
    row.names(arrays) <- 1:nrow(arrays)
    return(arrays)
    }
}



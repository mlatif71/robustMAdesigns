
  #
  GA.default <- function(dmat, layout, t, cmat, cinfo, leaveOut, popSize,
                         Pcross, Pmut, maxGen, verbose, crossType, selectType,
                         burnOut, convergeNum, nmarrays, keep.elite,
                         scaling, cname){
#######
   design <- dmat
   if (is.null(design)) {
      if (is.null(layout)) {
        message("design or layout must be specified!")
        return(NA)
      } else {
        design <- desMatrix(layout = layout)
      }
   }
   if (is.null(cinfo)) {
      if (is.null(dim(cmat)[1])) cinfo <- 1
      else cinfo <- dim(cmat)[1]
   }

   # convergence code : 0 -> perfect converging, 1 -> max and
   # median is too close degenerate pop, 2 -> reached max no of iteration
   converge <- 0
   keep.call <- match.call()
   design <- as.matrix(design)
   narray <- nrow(design)
   ncon <- length(cinfo)
   # generating initial population of designs of size ipopSize
   ii <- 1
   nFeas <-  0
   if (verbose) cat("\n")
   while (nFeas < floor(popSize*.5)) {
      ipop <- ini.pop(ii*popSize, design, t, cmat, cinfo, leaveOut,
                      nmarrays, scaling)
      nFeas <- nFeasible(ipop$opt)
      if (verbose) cat("   Number of feasible designs at initial population:",
                      nFeas, "\n")
      ii <- ii + 1
   }
   if (verbose) cat("\n")
   #############
   pop <-  ipop$pop
   opt <- ipop$opt
   opt.overall <- ipop$opt.overall
#   if(keep.elite) elite = list(pop=ipop$best.design, opt=ipop$max.opt,
#   opt.overall=ipop$max.opt.overall)
   #
   terminate <- 0
   trackMut <-  NULL
   des.h <- rbind(NULL, ipop$best.design)
   opt.h <- ipop$max.opt.overall
   if (keep.elite) elite <- list(design = ipop$best.design, opt = ipop$max.opt,
                       opt.overall = ipop$max.opt.overall)
   m <- stats::median(opt.overall, na.rm = T) # mean fitness value for the current generation
   s <-  stats::mad(opt.overall, na.rm = T) # sd of the fittness value for the current generation
   nFeas <-  nFeasible(opt)
   ###### print updates #############
   if (verbose)
      cat("   gen:  max(fitness)    median(fitness)    mad(fitness)   delta      nFeasible\n")
   if (verbose)
      cat("    ", 1,"        ", round(ipop$max.opt.overall, 4), "     ",
          round(m[1], 4), "     ", round(s[1], 3), "       ",
          round(abs(ipop$max.opt.overall - m[1]), 4), "     ",
          round(nFeasible(opt), 3), "\n")
   # generation loops
   elapsetime <- 0
   for (i in 2:maxGen) {
      t0 <- Sys.time()
      new.res <- evolution(pop, opt.overall, opt, cmat, cinfo, leaveOut,
                           design, Pcross, Pmut, popSize, crossType,
                           selectType, nmarrays, scaling)
      elapsetime[i] <-  difftime(Sys.time(), t0)
      ##########
      pop <-  new.res$pop
      opt <-  new.res$opt
      opt.overall <- new.res$opt.overall
      max.opt <-  new.res$max.opt
      max.opt.overall <-  new.res$max.opt.overall
      best.design <- new.res$best.design
      ####### keep elite ###############
      if (keep.elite) {
         if (new.res$max.opt.overall < elite$opt.overall) {
            max.opt.overall <- elite$opt.overall
            max.opt <- elite$opt
            best.design <- elite$design
            pop[popSize,] <- best.design
            opt[popSize,] <- max.opt
            opt.overall[popSize] <- max.opt.overall
            elite <- list(design = best.design, opt = max.opt,
                          opt.overall = max.opt.overall)
         }
         else {
            elite <- list(design = best.design, opt = max.opt,
                          opt.overall = max.opt.overall)
         }
      } # keep.elite
      #################
      des.h <- rbind(des.h, best.design)
      opt.h <- c(opt.h, max.opt.overall)
      m <- c(m, stats::median(opt.overall, na.rm = T))
      s <- c(s, stats::mad(opt.overall, na.rm = T))
      nFeas <- c(nFeas, nFeasible(opt))
      #
      ######### convergence criterion ############
      if (i >= burnOut) {
         ind <- opt.h[i] - opt.h[i - 1]
         if (abs(ind) < .001) terminate <- terminate + 1
         else terminate <- 0
      }
      if (abs(opt.h[i] - m[i]) < .00001) terminate <- 100
      if (s[i] <= .0001 | nFeas[i] < popSize / 10) terminate <- 100
      converge <- 1
      ######## if convergence criterion satisfies #################
      if (terminate >= convergeNum) {
         if (terminate == 100) {
            if (verbose) cat("   Process terminates because maximum and average of fitness value are too close! \n\n")
         }
         if (verbose) cat("\n")
         if (ncon > 1) opt <- round(opt, 3)
         else opt <- round(opt.overall, 3)
         opt.h <- round(opt.h, 3)
         best.design <- des.h[which.max(opt.h), ]
         opt.best.design <- eCriteria(dmat = design[best.design, ],
                                      cmat = cmat,
                                      cinfo = cinfo, type = "e")$eff
         ######### output ######
         #
         ## bring opt.overall into its own scale
         opt.overall <- apply(opt, 1, stats::weighted.mean, leaveOut,
                            na.rm = T)
         #
         history <- list(design = des.h, opt = opt.h, m = m, s = s,
                         nFeas = nFeas, nGen = i, Pmut = eval(Pmut),
                         Pcross = eval(Pcross), elapse.time = elapsetime)
         ##
         best <- list(design = best.design, opt = opt.best.design,
                      opt.overall = max(opt.h, na.rm = T))
         input <- list(dmat = design, cmat = cmat, cinfo = cinfo,
                       call = keep.call)
         #
         if (!is.null(cname)) colnames(opt) <- cname
         out <-  list(design = pop, opt.overall = opt.overall, opt = opt,
                      best = best,
                       history = history, converge = converge, input = input)
         return(out)
      } # terminate
      ##########################
      if (verbose)
        cat("    ", i,"        ", round(max.opt.overall, 4), "     ",
            round(m[i], 4), "     ", round(s[i], 3), "       ",
            round(abs(max.opt.overall - m[i]),4), "      ",
            round(nFeas[i], 3), "\n")
   } # maxGen
   if (verbose) cat("\n")
   if (verbose) cat("Maximum number of generations achieved!\n")
   best.design <- des.h[which.max(opt.h), ]
   opt.best.design <- eCriteria(dmat = design[best.design,],
                                cmat = cmat,
                                cinfo = cinfo, type = "e")$eff
   best <- list(design = best.design, opt = opt.best.design,
                opt.overall = max(opt.h, na.rm = T))
   history <- list(design = des.h, opt = opt.h, m = m, s = s,
                   nFeas = nFeas, nGen = maxGen)
   input <-  list(dmat = design, cmat = cmat, cinfo = cinfo,
                  call = keep.call)
   opt.overall <-  apply(opt, 1, stats::weighted.mean, leaveOut, na.rm = T)
   if (!is.null(cname)) colnames(opt) <- cname
   out <-  list(design = ipop$pop, opt.overall = opt.overall, opt = opt,
               best = best, history = history, converge = 2, input = input)
   return(out)
} # GA

# estimating optimality values for given set of designs
# if no designs are given it generate a set of designs
fitness <-  function(ipop, cmat, cinfo, leaveOut, design, nmarrays, scaling){
   #
   # number of contrast
   ncon <- length(cinfo)
   #
   # compute optimality values for cmat and design
   opt <- t(apply(ipop, 1,
                  function(x) rfitness(
                    dmat = design[x,],
                    cmat = cmat,
                    cinfo = cinfo,
                    type = "e",
                    nmarrays = nmarrays)))
   colnames(opt) <- NULL
   #
   # mean of the optimality value over all contrasts
   if (ncon > 1) {
      opt.overall <- apply(opt, 1, stats::weighted.mean, leaveOut, na.rm = T)
      #
      # if all effects are nonestimable then the weighted  mean produce nan
      opt.overall[is.nan(opt.overall)] <- NA
      #
      ##### scaling ###### Ref. Greenberg (89), p. 76
#     if(scaling){
#        avg = mean(opt.overall, na.rm=T)
#        mx = max(opt.overall, na.rm=T)
#        mn = min(opt.overall, na.rm=T)
#        if(mn > 2*avg - mx){
#           delta = mx - avg
#           b = avg*(mx - 2*avg)/delta
#        }
#        else {
#           delta = avg - mn
#           b = -mn*avg/delta
#        }
#        a = avg/delta
#        opt.overall = a*opt.overall + b
#     } # scaling
      #
      opt.overall <- round(opt.overall, 6)
      #
      # ordering wrt opt.overall
      #
      ipop <- ipop[order(opt.overall, decreasing = T),]
      opt <- opt[order(opt.overall, decreasing = T), ]
      opt.overall <- sort(opt.overall, decreasing = T, na.last = TRUE)
      ##
      # find designs that can estimalte all the contrasts
      indx <- !is.na(apply(opt, 1, sum))
      sindx <- sum(indx)
      #
      ########## best one ########
      if (sindx == 0) {
        max.opt.overall <- 0
        max.opt <- rep(0, ncon)
        best.design <- rep(0, dim(design)[2])
      }
      else {
        opt1 <- opt[indx,]
        opt.ov1 <- opt.overall[indx]
        ipop1 <- ipop[indx,]
        mx <- which.max(opt.ov1)[1]
        ##
        if (sindx == 1) {
          max.opt <- opt[indx,]
          best.design <- ipop[indx,]
        }
        else {
          max.opt <- opt1[mx,]
          best.design <- ipop1[mx,]
        }
        max.opt.overall <- stats::weighted.mean(max.opt, leaveOut, na.rm = T)
     }
   } # ncon > 1
   else {
      opt <- c(opt)
      opt[is.na(opt)] <- 0
      ipop <- ipop[order(opt, decreasing = T), ]
      opt.overall <- sort(opt, decreasing = T)
      opt <- opt.overall
      #### best one #########
      mx <- which.max(opt)
      max.opt <- opt[mx]
      max.opt.overall <- max.opt
      best.design <- ipop[mx, ]
   } # else ncon
   #
   return(list(opt = opt, opt.overall = opt.overall, pop = ipop,
               max.opt = max.opt,
                 max.opt.overall = max.opt.overall, best.design = best.design))
}

# generating initial population
ini.pop <- function(popSize, design, t, cmat, cinfo, leaveOut, nmarrays, scaling){
   nr <- nrow(design)
   ncon <- length(cinfo)
   # genereting initial designs
   ipop <- matrix(sample(1:nr, size = t * popSize, replace = T),
                  ncol = t, byrow = TRUE)
   ####
   out <- fitness(ipop, cmat, cinfo, leaveOut, design, nmarrays, scaling)
   return(list(pop = out$pop, opt = out$opt, opt.overall = out$opt.overall,
               best.design = out$best.design,
                      max.opt = out$max.opt,
               max.opt.overall = out$max.opt.overall))
}

# evolution
evolution <- function(pop, opt, optAll, cmat, cinfo, leaveOut, design,
                      Pcross, Pmut,
                        popSize, crossType, selectType, nmarrays, scaling){
   t <- ncol(pop)
   ncon <- length(cinfo)
  # s1 = summary(opt)
   if ( ncon > 1) {
      #ind = apply(!is.na(optAll), 1, sum)
      ind <- apply(is.na(optAll), 1, sum)
      tem <- tapply(opt, ind, min)
      ltem <- length(tem)
      if (ltem > 1) {
        for (ii in 2:ltem) {
           nn <- as.numeric(names(tem)[ii])
           nn0 <- as.numeric(names(tem)[ii - 1])
           opt[ind == nn & opt > tem[ii - 1]] <- max(
             tem[ii - 1] - (nn0 + 1) * tem[ii - 1] / ncon, 0)
            #+ opt[ind==nn & opt>tem[ii-1]]/100
           tem[ii] <- min(opt[ind == nn])
        }
      }
   }
#   else ind = as.numeric(!is.na(optAll))
   #s2 = summary(opt)
   # penalty function
#   opt = (ind/ncon)*opt
#
####
   if (scaling) opt <- Scaling(opt)
 #  s3  = summary(opt)
 #  print(rbind(s1,s2,s3))
##
   opt <- round(opt, 6)
   spop1 <- recombination(pop, opt, nrow(design), Pcross, Pmut, crossType, selectType, popSize)
   spop <- spop1$pop
   out <- fitness(ipop = spop, cmat, cinfo, leaveOut, design, nmarrays, scaling)
   return(list(pop = out$pop, opt = out$opt, opt.overall = out$opt.overall,
               nMut = spop1$nMut, nCross = spop1$nCross,
               Pmut = Pmut, best.design = out$best.design,
               max.opt = out$max.opt,
               max.opt.overall = out$max.opt.overall))
}

# recombination step contains (i) crossover and (ii) mutation
recombination <- function(pop, opt, narrays, Pcross, Pmut, crossType,
                          selectType, popSize){
# input
#       pop : population of current generation
#       narrays : number of possible arrays
#       sprob : selection prob for each individual
#       Pcross : crossover prob
#       Pmut : mutation prob
#
   t <- ncol(pop)
   Size <- length(opt)
   if (selectType == "SPF") {
      # first select popsize individual for next generation wrt the selection prob
      opt1 <- opt
      opt1[is.na(opt)] <- 0
      #min(opt1, na.rm=T)/2
      tem <- sample(1:Size, size = popSize, replace = T, prob = opt1)
      # pair the selected individuals for mating
      spair <- matrix(sample(tem, size = popSize, replace = F),
                      ncol = 2, byrow = T)
   }
   else if (selectType == "RSS") {
      rssRes <- RemStochSamp(pop = pop, opt = opt, popSize)
      pop <- rssRes$pop
      spair <- matrix(sample(1:nrow(pop), size = popSize, replace = F),
                      ncol = 2, byrow = T)
   }
   else{
      stop("This method of selection operator has not been implemented yet.\n")
   }
   # decide which pair perform the crossover
   t1 <- Pcross
   t0 <- (1 - Pcross)
   pcross <- sample(x = c(1, 0), size = popSize / 2,
                    replace = T, prob = c(t1, t0))
   #number of crossover
   nCross <- sum(pcross)
   ### mutation
   pmut <- sample(x = c(1,0), size = length(pop),
                  replace = T, prob = c(Pmut, 1 - Pmut))
   #pmut = sample(x=c(1,0), size=popSize*t, replace=T, prob=c(Pmut, 1-Pmut))
   pmut <- matrix(pmut, ncol = t, byrow = T)
   res <- lapply(1:nrow(spair), function(i){
             if (pcross[i]) { # crossover
                return(crossover(pop[spair[i, ], ],
                                 crossType,
                                 pmut[spair[i, ], ],
                                 narrays, t))
             }# no crossover
             else (return(pop[spair[i, ], ]))
             #else (return(list(out=pop[spair[i,],], nMut=0)))
          })
   #
    out <- as.matrix(do.call("rbind", res))
   row.names(out) <- NULL
   return(list(pop = out, nMut = sum(pmut), nCross = nCross))
}
#
# crossover operator
crossover <- function(d, type, mutation = NULL, narrays = NULL, ld){
#   ld = ncol(d)
   d0 <- d
   if (type == "one.point") {
#      print(d0)
      xx <- sample(1:(ld - 1), size = 1)
      d0[1, ] <- c(d[1, 1:xx], d[2, (xx + 1):ld])
      d0[2, ] <- c(d[2, 1:xx], d[1, (xx + 1):ld])
     # return(d0)
   }
   else if (type == "two.points") {
      xx <- sort(sample(1:(ld - 1), size = 2, replace = F))
      d0[1, ] <- c(d[1, 1:xx[1]], d[2, (xx[1] + 1):xx[2]],
                   d[1, (xx[2] + 1):ld])
      d0[2,] <- c(d[2,1:xx[1]], d[1, (xx[1] + 1):xx[2]],
                  d[2, (xx[2] + 1):ld])
     # return(d0)
   }
   else if (type == "uniform") {
      xx <- matrix(sample(c(0,1), size = 2*ld, replace = T),
                   ncol = ld, byrow = T)
      for (i in 1:2) {
        d0[i,which(xx[i,] == 1)] <- d[1, which(xx[i,] == 1)]
        d0[i,which(xx[i,] == 0)] <- d[2, which(xx[i,] == 0)]
      }
      #return(d0)
   }
   else stop("This method has not been implemented yet\n")
   # mutation step
      #print(d0)
   if (is.null(mutation) | sum(mutation) == 0) return(d0)
   else {
      d1 <- sapply(1:2,
              function(i){
                aArrays <- (1:narrays)[-d0[i,]]
                if (length(aArrays) == 0) aArrays <- 1:narrays
                sd0 <- sum(mutation[i,])
                if (sd0 == 0) return(d0[i,])
                else {
                  rArrays <- sample(aArrays, size = sd0, replace = T)
                  d0[i, which(mutation[i,] > 0)] <- rArrays
                  return(d0[i,])
                }
          })
   }
   return(t(d1))
}

# mutation
mutation <- function(pop, design){
   ind1 <- sample(1:length(pop), size = 1, replace = F)
   ind2 <- sample(1:nrow(design), size = 1, replace = F)
   pop[ind1] <- ind2
   return(pop)
}

# Remainder Stochastic Sampling
# Procedure :
#             (i) calculate the number of expected frequencies (fi/\bar f) from the given optimality values
#             (ii) integer part of the expected frequencies indicate the number of the times the individual is copied in the next
#             generation
#             (iii) the decimal part of the expected frequencies is used as prob to sample the remaining number of the individuals in
#             the next generation
#
RemStochSamp <- function(opt, pop, nsize){
   nr <- nrow(pop)
   nexf <- length(opt)*(opt/sum(opt, na.rm = T))
   #
   #### what to do with designs whose expected freq is NA #######
   # assign a very small value
   #
   nexf[is.na(nexf)] <- 0
   #min(nexf, na.rm=T)/2
   dup <- floor(nexf)
   if (length(dup) >= 1) {
     probs <- nexf - dup
     tem1 <- rep(1:nr, dup)
     remNo <- nsize - length(tem1)
   }
   else {
     remNo <- nsize
     tem1 <- NULL
   }
   if (remNo > 0) {
     tem2 <- sample(1:nr, size = remNo, prob = probs)
   }
   else tem2 <- NULL
   tem <- c(tem1, tem2)
   npop <- pop[tem,]
   nopt <- opt[tem]

   return(list(pop = npop))
}

### Number of feasible designs
nFeasible = function(opt){
   if (is.matrix(opt)) return(sum(!is.na(apply(opt, 1, sum))))
   else return(sum(!is.na(opt)))
}

########### average efficiency in the presence of missing arrays

# this codes computes robustness criteria
# Mahbub Latif: Novermber 08, 2004
rfitness <- function(dmat, cmat, cinfo=NULL, nmarrays=0, type="e"){
  #
   if ((dim(dmat)[1] - nmarrays) < 2) {
     message("For this design matrix \"nmarrays\" must be < ",
          dim(dmat)[1] - 1)
     return(NA)
   }
  #
   if (is.null(cinfo)) {
      if (is.null(dim(cmat))) cinfo <- 1
      else cinfo <- dim(cmat)[1]
   }
  #
   nc <-  length(cinfo)
   frun <-  unlist(eCriteria(dmat = dmat, cmat = cmat, m = 0,
               cinfo = cinfo, type = type)$eff)
   if (nmarrays == 0) return(frun)
   pED <- rbind(NULL, !is.na(frun))
   avgEff <-  frun
   ex.sl <- 1
   while (ex.sl <= nmarrays) {
      res <- unlist(eCriteria(dmat = dmat, cmat = cmat, m = ex.sl,
                cinfo = cinfo, type = type)$eff)
      if (nc > 1) nr.res <- nrow(res)
      else nr.res <- length(res)
#
      if (nc > 1) {
          pED  <-  apply(res, 2, function(x) sum(!is.na(x)))/nr.res
          if (sum(pED) >= nc) avgEff <- avgEff +
              apply(res, 2, sum, na.rm = TRUE) / nr.res
      }
      else {
         pED <- sum(!is.na(res))/nr.res
         if (pED >= 1) avgEff <- avgEff + sum(res, na.rm = TRUE)/nr.res
      }
#
      ex.sl <- ex.sl + 1
   }
#
   names(avgEff) <- NULL
   return(avgEff/(nmarrays + 1))
}
#######
#
# scaling
Scaling <- function(opt.overall){
   avg <- mean(opt.overall, na.rm = T)
   mx <- max(opt.overall, na.rm = T)
   mn <- min(opt.overall, na.rm = T)
   if (mn > 2*avg - mx) {
      delta <- mx - avg
      b <- avg*(mx - 2*avg)/delta
   }
   else {
      delta <- avg - mn
      b <- -mn*avg/delta
   }
   a <- avg/delta
   opt.overall <- a*opt.overall + b
   return(opt.overall)
}
#####

Fun1.eCriteria <- function(x=NULL, design.list, cmat, cinfo, tol, type){
    cinfo <- cumsum(cinfo)

    if (is.null(x)) design.list.1 <- design.list
    else{
       design.list.1 <- design.list
       design.list.1 <- design.list.1[-c(x), ]
    }
    res <- NULL
    k <- 1
    for (i in cinfo) {
       if (length(cinfo) == 1) tem <- cmat # 02.02.04
       else {
         tem <- cmat[k:i, ]
         }
       if (is.null(dim(tem))) tem <- t(as.matrix(tem))
       res0 <- varFac(x = design.list.1, cc = tem,
                      tol = tol, type = type)
       res <- rbind(res, res0)
       k <- i + 1
    }
   if (is.null(x)) return(res)
   else return(c(res, d1 = x))
}


#estimable <- function(dmat, cmat) {
#   res <- varFac(dmat, cmat, 1.e-5, "e")
#   if (is.na(res)) return(FALSE)
#   else return(TRUE)
#}


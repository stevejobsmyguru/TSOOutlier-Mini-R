
##FIXME
# see add argument "duplicates.remove = TRUE" as an option 
# for non-automatic inspection of the procedure

# types of outliers
# types = c("IO", "AO", "LS", "TC", "SLS")

locate.outliers <- function(resid, pars, cval = 3.5, 
  types = c("AO", "LS", "TC"), delta = 0.7, n.start = 50)
{
  # mean absolute deviation of residuals

  sigma <- 1.483 * quantile(abs(resid - quantile(resid, probs = 0.5, na.rm = TRUE)), 
    probs = 0.5, na.rm = TRUE)
  #n <- length(resid)

  # t-statistics (and coefficents) at every time point for the 
  # selected types of outliers

  tmp <- outliers.tstatistics(pars = pars, resid = resid, 
    types = types, sigma = sigma, delta = delta, n.start = n.start)

  #ind <- which(abs(tmp[,,"tstat"]) > cval, arr.ind = TRUE)
  #mo <- data.frame(
  #  factor(gsub("^(.*)tstats$", "\\1", dimnames(tmp)[[2]][ind[,2]]),
  #    levels = c("IO", "AO", "LS", "TC", "SLS")), ind[,1],
  #  tmp[,,"coefhat"][ind], tmp[,,"tstat"][ind])
  #colnames(mo) <- c("type", "ind", "coefhat", "tstat")
  #if (nrow(ind) == 1) # otherwise the row name is "row" instead of "1"
  #  rownames(mo) <- NULL
  ind <- which(abs(tmp[,,"tstat",drop=FALSE]) > cval, arr.ind = TRUE)
  mo <- data.frame(
    factor(gsub("^(.*)tstats$", "\\1", dimnames(tmp)[[2]][ind[,2]]),
      levels = c("IO", "AO", "LS", "TC", "SLS")), ind[,1],
    tmp[,,"coefhat",drop=FALSE][ind], tmp[,,"tstat",drop=FALSE][ind])

  colnames(mo) <- c("type", "ind", "coefhat", "tstat")
  if (nrow(ind) == 1) # otherwise the row name is "row" instead of "1"
    rownames(mo) <- NULL

  #changed in version 0.6-4: 
  #the following is moved to locate.outliers.iloop(), find.consecutive.outliers()
  #
  # remove consecutive LS
  #
  # if a LS is found at consecutive time points, then 
  # keep only the point with the highest abs(tstat)
  #
  #resLS <- mo[mo[,"type"] == "LS",] 
  #nls <- nrow(resLS)    
  #if (nls > 1)  # && any(diff(sort(resLS[,"ind"])) == 1)
  #{
  #  #changed in version 0.6-4
  #  resLS <- resLS[order(resLS[,"ind"]),]
  #  ids <- c(0, which(diff(resLS[,"ind"]) != 1), nls)
  #  #apparently here 'lapply' is slightly faster than 'sapply'
  #  aux <- lapply(as.list(seq(along = ids[-1])), 
  #    function(i, x) x[seq.int(ids[i]+1, ids[i+1])], x = seq_len(nls))
  #  rmid <- NULL
  #  for (i in seq_along(aux))
  #  {
  #    if (length(aux[[i]]) == 1)
  #      next 
  #    id <- which.max(abs(resLS[aux[[i]],"tstat"]))
  #    id <- aux[[i]][-id]
# #debug
#stopifnot(length(id) > 0)
  #    rmid <- c(rmid, id)
  #  }
  #  if (length(rmid) > 0)
  #  {    
  #    #changed in version 0.6-4 
  #    #mo <- mo[-as.numeric(rownames(resLS[rmid,])),]
  #    mo <- mo[-match(as.numeric(rownames(resLS[rmid,])), rownames(mo)),]
  #  }
  #}

  # deal with multiple types of outliers at the same time point (if any)
  # if the t-statistics of more than one type outlier exceed "cval"
  # at a given time point, then keep the outlier with the highest abs(tstat)

#changed in version 0.6-5: added 'unique'
  ref <- unique(mo[,"ind"][duplicated(mo[,"ind"])])
  for (i in ref)
  {
    ind <- which(mo[,"ind"] == i)
    moind <- mo[ind,]

    ##NOTE
    # in a previous version precedence was given to AO, LS and TC over IO; 
    # in that case, an IO was found only when the t-statistic of the IO 
    # exceeded the threshold "cval" and when the t-statistics for the other 
    # types of outliers were below the critical value;
    # if that rule is applied again in a future version, check that 
    # this approach to remove duplicates is compatible with it

    #if ("IO" %in% moind[,"type"]) {
    #  tmp <- moind[which(moind[,"type"] != "IO"),]
    #} else
    #  tmp <- moind[which.max(abs(moind[,"tstat"])),]

#changed in version 0.6-5 (simplified)
    #tmp <- moind[which.max(abs(moind[,"tstat"])),]
    #mo <- mo[-ind,]
    #mo <- rbind(mo, tmp)
    mo <- mo[-ind[-which.max(abs(moind[,"tstat"]))],]
  }

  mo
}

locate.outliers.iloop <- function(resid, pars, cval = 3.5, 
  types = c("AO", "LS", "TC"), maxit = 4, delta = 0.7, n.start = 50,
  logfile = NULL)
{
  if(!is.ts(resid))
    stop(paste(sQuote("resid"), "must be a", sQuote("ts"), "object"))

  n <- length(resid)
  s <- frequency(resid)
##FIXME see comment this line out
  moall <- data.matrix(numeric(0)) # nrow(moall) == 0 #NULL
  iter <- 0

  # begin inner loop

  while (iter < maxit)
  {
    mo <- locate.outliers(resid = resid, pars = pars, cval = cval, 
      types = types, delta = delta, n.start = n.start)

    if (!is.null(logfile))
    {
      msg <- paste("\niloop, iteration:", iter, "\n")
      cat(msg, file = logfile, append = TRUE)
      capture.output(mo, file = logfile, append = TRUE)
    }

    cond <- nrow(mo) > 0

    # discard outliers identified at consecutive time points (if any)
    #
    #changed in version 0.6-5: 
    # this cleaning is now done here and in locate.outliers.oloop() instead 
    # of in locate.outliers()
    # in addition to LS, other types of outliers are also checked 

    if (cond)
    {
      rmid <- c(
        find.consecutive.outliers(mo, "IO"),
        find.consecutive.outliers(mo, "AO"),
        find.consecutive.outliers(mo, "LS"),
        find.consecutive.outliers(mo, "TC"),
        find.consecutive.outliers(mo, "SLS"))

      if (length(rmid) > 0)
      # do not use is.null(rmid), 'rmid' may be NULL or 'character(0)'
      {
         mo <- mo[-rmid,]
      }
    }

    # remove duplicates in 'moall' and in 'mo' (last and current iterations)
    #
    # remove outliers at those points where an outlier was found in a previous 
    # iteration within this loop;
    # for example, if an AO was found at observation 15, it was removed from
    # the residuals according to "outliers.regressors", then we do not expect the 
    # detection of another AO or other type of outlier at the same time point;
    # if that happens, the type of outlier detected the first time is kept;
    # example: duplicates at this point are found in series log(hicp[["000000"]])
    # in that case this loop takes 2 iterations when duplicates are removed, 
    # without removing duplicates it takes 4 iterations

    #if (remove.duplicates)
    if (cond && iter > 0)
    {
      id.dups <- na.omit(match(moall[,"ind"], mo[,"ind"]))
      if (length(id.dups) > 0)
      {
        mo <- mo[-id.dups,]
        cond <- nrow(mo) > 0
      }
      # no problems with mo[-id.dups,]
      # the two dimensions of a matrix are kept in data.frame "mo" even if 
      # "mo" contains one element after removing duplicates
    }

    if (!cond)
      break

    moall <- rbind(moall, mo)

    ##NOTE
    # at this point, outliers at consecutive time points may arise in 'moall', 
    # however, it may be better not to deal with here;
    # if after checking their t-statistics it is concluded that an outlier 
    # from a previous iteration should now be removed, it won't be easy 
    # to restore the vector 'resid' that was modified below 
    # (the outliers and 'pars0 from each iteration should be stored); 
    # it seems simpler to deal with this in locate.outliers.oloop()
    #
    # this is not an issue with possible duplicates in 'moall' and 'mo';
    # as done above, the outlier from the previous iteration is kept, hence, 
    # only the current 'mo' would be modified before updating the vector 'resid'

    oxreg <- outliers.regressors(pars = pars, mo = mo, n = n, weights = TRUE,
      delta = delta, freq = s, n.start = n.start)

    #resid0 <- resid
    resid <- resid - rowSums(oxreg)

    iter <- iter + 1
  } # end while

  # in practice this warning may be ignored or inspected further 
  # for example incrementing the value of "maxit"

  if (iter == maxit)
    warning(paste("stopped when", sQuote("maxit.iloop"), "was reached"))

##FIXME
# see return "iter"

  moall
}

locate.outliers.oloop <- function(y, fit, types = c("AO", "LS", "TC"), 
  cval = NULL, maxit.iloop = 4, maxit.oloop = 4, delta = 0.7, n.start = 50, logfile = NULL)
{
  n <- length(y)
  s <- frequency(y)

  # when "fit" is the output of "auto.arima" the argument "y" 
  # could be avoided and set "y <- fit$x"
  # but if "fit" is created by "arima" there is no way to 
  # get the series "y"

  # default critical value 
  # (same as in functions "tso" and "locate.outliers.oloop")

  if (is.null(cval))
  {
    if (n <= 50) {
      cval <- 3
    } else 
    if (n >= 450) {
      cval <- 4
    } else
      cval <- round(3 + 0.0025 * (n - 50), 2)
  }

  # tail(): take the last element just in case fit$call[[1]] is 
  # for example "forecast::auto.arima"
  tsmethod <- ifelse(inherits(fit, "stsmFit"), 
    "stsm", tail(as.character(fit$call[[1]]), 1))
  #s <- frequency(y)
  moall <- data.frame(matrix(nrow = 0, ncol=4, 
    dimnames = list(NULL, c("type", "ind", "coefhat", "tstat"))))
  iter <- 0

  # index of initial residuals

  if (inherits(fit, "Arima")) {
    tmp <- fit$arma[6] + fit$arma[5] * fit$arma[7]
    id0resid <- if (tmp > 1) seq.int(tmp) else c(1, 2)
  } else 
  if (inherits(fit, "stsmFit")) {
    id0resid <- seq_len(n - length(fit$model@diffy))
  } else
    stop("unexpected type of fitted model")

  # begin outer loop

  #while (TRUE)
  while (iter < maxit.oloop)
  {
    # extract the necessary information from the fitted model
    # parameter estimates and residuals

    pars <- switch(tsmethod, 
      "auto.arima" = , "arima" = coefs2poly(fit),
      #"stsm" = stsm::char2numeric(fit$model)) #Hari is trying undok stsm to
      "stsm" = char2numeric(fit$model)) #Hari is trying undok stsm to

    ##NOTE by default residuals(fit, standardised = FALSE) 
    # only relevant for "stsm" but the argument could set here 
    # explicitly, it would be ignored if "fit" is an "Arima" object
    resid <- residuals(fit)

    # impute NAs in order to avoid propagation of NAs 
    # in the t-statistics computed by 'outliers.tstatistics'
    idna <- which(is.na(resid))
    if (length(idna))
      resid[idna] <- mean(resid, na.rm=TRUE)
  
    if (any(abs(na.omit(resid[id0resid])) > 3.5 * sd(resid[-id0resid], na.rm = TRUE)))
    {
##FIXME
# see add factor 3.5 as argument

      # this was necessary since in one series (I think it was hicp[["000000"]])
      # the first residuals were too erratic and caused the procedure to detect
      # too many outliers in the whole series

      resid[id0resid] <- 0
      warning(paste("the first", tail(id0resid, 1), "residuals were set to zero"))
      
      ##NOTE
      # if this warning is returned
      # the first observations of the series may need to be inspected
      # for possible outliers, since those observations were ignored by the procedure
    }

    # locate possible outliers

    mo <- locate.outliers.iloop(resid = resid, pars = pars, cval = cval, 
      types = types, maxit = maxit.iloop, delta = delta, n.start = n.start, 
      logfile = logfile)

    # discard outliers identified at consecutive time points (if any)
    #
    # this is done by type of outlier, e.g., two or more consecutive LS are 
    # replaced by one LS (the one with the highest abs(t-statistic);
    # it is still possible to get two or more consecutive outliers of different 
    # type at consecutive time points; 
    # alternatively consecutive outliers of any type could be replaced by a 
    # single outlier. For the time being, this is kept in this way so that 
    # for example, the following sequence AO1,LS2,LS3,LS4 can collapse to 
    # AO1,LS4, i.e., the AO is kept
    #
    # it is not enough to do this in locate.outliers.iloop(), because 
    # only the outliers detected at a given interation are checked there;
    # here, the whole set of outliers returned by locate.outliers.iloop() 
    # is checked

    if (nrow(mo) > 0)
    {
      rmid <- c(
        find.consecutive.outliers(mo, "IO"),
        find.consecutive.outliers(mo, "AO"),
        find.consecutive.outliers(mo, "LS"),
        find.consecutive.outliers(mo, "TC"),
        find.consecutive.outliers(mo, "SLS"))
      #rmid <- NULL
      #for (type in types)
      #  rmid <- c(rmid, find.consecutive.outliers(mo, type))

      if (length(rmid) > 0)
      # do not use is.null(rmid), 'rmid' may be NULL or 'character(0)'
      {
        #changed in version 0.6-4
        #mo <- mo[-as.numeric(rownames(mo[rmid,])),]
        #mo <- mo[-match(as.numeric(rownames(mo[rmid,])), rownames(mo)),]
        #changed in version 0.6-5
        mo <- mo[-rmid,]
      }
    }

    # remove duplicates (if any)
    # similar to locate.outliers.iloop(), if an outlier is detected at 
    # an observation where some type of outliers was already detected in 
    # a previous run, the outlier that was detected first is kept

    if (nrow(mo) > 0 && iter > 0)
    {
      id.dups <- na.omit(match(moall[,"ind"], mo[,"ind"]))
      if (length(id.dups) > 0)
        mo <- mo[-id.dups,]
      # no problems with mo[-id.dups,]
      # the two dimensions of a matrix are kept in data.frame "mo" even if 
      # "mo" contains one element after removing duplicates
    }

##FIXME
#at this point it could be checked for consecutive outliers of a given type
#and keep the one detected in a previous iteration
#if (nrow(mo) > 0 && iter > 0)
#tmp <- rbind(mo, moall)
#rmid <- c(
#  find.consecutive.outliers(tmp, "IO"),
#  find.consecutive.outliers(tmp, "AO"),
#  find.consecutive.outliers(tmp, "LS"),
#  find.consecutive.outliers(tmp, "TC"),
#  find.consecutive.outliers(tmp, "SLS"))
#print(rmid)

    if (!is.null(logfile))
    {
      msg <- paste("\noloop, iteration:", iter, "\n")
      cat(msg, file = logfile, append = TRUE)
      capture.output(mo, file = logfile, append = TRUE)
    }

    # this must be here, not before the previous "if" statement
    # since it may modify "mo"

    if (nrow(mo) == 0)
      break

    moall <- rbind(moall, mo)

    # remove the effect of outliers on the data and 
    # fit the model for the adjusted series

    oeff <- outliers.effects(mo = mo, n = n, weights = TRUE, 
      delta = delta, pars = pars, n.start = n.start, freq = s)

    # 'y' is overwritten; 'oeff' is based on 'mo' not 'moall'

    y <- y - rowSums(oeff)

    switch(tsmethod,
##FIXME 
#if 'fit' includes intercept or drift, pass here (it could be done based on names of coef(fit))

       # do not modify and evaluate the call, i.e. do not run eval(fit$call)
       # since it will run the model selection procedure (if tsmethod = "auto.arima")
       # here we only want to refit the model (not choose or select a model)

       "auto.arima" = fit <- arima(y, order = fit$arma[c(1,6,2)], 
         seasonal = list(order = fit$arma[c(3,7,4)])),

      # this reuses arguments passed to the optimization method, e.g. method = "CSS",
      # if they were specified when the input object "fit" was created,
      # (for example through argument "args.tsmethod" of function "tso")

      "arima" = {
        fitcall <- fit$call
        fitcall$x <- y
        # rename fitcall$series since it can be a very long character string to store
        #fitcall$series <- "x"
        fit <- eval(fitcall)
      },

      "stsm" = {
        fitcall <- fit$call
        ##NOTE
        # fitcall$x contains the model, not fitcall$m, since now "stsmFit" is called 
        # instead of "maxlik.td.optim" and the other functions
        fitcall$x@y <- y
        dy <- fitcall$x@fdiff(y, frequency(y))
        fitcall$x@diffy <- dy
        if (!is.null(fitcall$x@ssd))
          fitcall$x@ssd <- Mod(fft(as.numeric(dy)))^2 / (2*pi*length(dy))
        ##NOTE
        #last parameter estimates, fit$pars, could be used as starting values
        #fitcall$x@pars[] <- fit$pars
        fit <- eval(fitcall)
      }
    )

    if (!is.null(logfile))
    {
      msg <- paste("\nmodel chosen and fitted for the adjusted series:\n")
      cat(msg, file = logfile, append = TRUE)
      capture.output(fit, file = logfile, append = TRUE)
    }

    iter <- iter + 1
  } # end while

  if (iter == maxit.oloop)
    warning(sprintf("stopped when \'maxit.oloop = %d\' was reached", maxit.oloop))

  # time points with multiple types of potential outliers are not expected
  # since they are removed at each iteration
  # stopifnot(!any(duplicated(moall[,"ind"])))  

  if (any(duplicated(moall[,"ind"])))
  {
    # stop for debugging 
    # so far this event has not occurred
    stop("unexpected duplicates since they are handled within the loop above")
  }

  # "coefs" is not actually used but keep so far, it may used to see 
  # the estimates for external regressors "xreg" which are not included in "pars"

  list(fit = list(coefs = coef(fit), pars = pars,
    resid = resid, n = n), outliers = moall, iter = iter)
}

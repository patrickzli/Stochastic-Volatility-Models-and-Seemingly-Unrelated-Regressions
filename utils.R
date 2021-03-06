# utility functions that will be helpful for this course

#' Numerical check for MLE.
#'
#' Given a log-likelihood function and a potential MLE, checks whether each one-dimensional version of the log-likelihood which varies one parameter at a time with all others at the MLE is indeed maximized at the MLE value.
#'
#' @param loglik loglikelihood function.  Takes a single vector argument.
#' @param theta.mle The potential MLE.
#' @param itheta indices of one dimensional functions to evaluate and plot.  Defaults to all parameters.
#' @param theta.names Optional vector of parameter names for plotting.
#' @param theta.rng Optional two-row matrix giving the plotting limits for each parameter.  Defaults to theta.mle +/- .5 * abs(theta.mle)
#' @param refit If \code{TRUE}, recalculates the range so that drop is more or less the same on either side of \code{theta.mle}.
#' @param layout Optional vector giving the number of rows and columns in the plot.  For \code{p} parameters, defaults to \code{c(nr, nc)}, where \code{nr = floor(p)} and \code{nc = ceiling(p/nr)}.
#' @param plot Logical, whether or not to plot the conditional likelihoods.  See Value for details.
#' @return A logical vector of length \code{p}, indicating whether or not the MLE was the maximum of the plotted points along each coordinate.  This is for numerically checking what the plots otherwise do visually.
mle.check <- function(loglik, theta.mle, itheta, theta.names, theta.rng,
                      refit = TRUE, layout, plot = TRUE) {
  ntheta <- length(theta.mle) # number of parameters
  npts <- 200 # number of points to plot
  is.mle <- rep(NA, ntheta) # check whether MLE maximizes 1d loglikelihoods
  ll.max <- loglik(theta.mle) # maximum value
  if(missing(itheta)) itheta <- 1:ntheta
  if(is.logical(itheta)) itheta <- which(itheta) # convert T/F's to indices
  if(missing(theta.names)) {
    theta.names <- paste0("theta[",1:ntheta,"]")
    # converts to expression so symbol "theta_i" is plotted
    theta.names <- parse(text = theta.names)
  }
  if(missing(theta.rng)) {
    theta.rng <- rbind(theta.mle - .5 * abs(theta.mle),
                       theta.mle + .5 * abs(theta.mle))
  }
  # shorten theta.names and theta.rng if necessary
  ntheta2 <- length(itheta)
  if(length(theta.names) > ntheta2) theta.names <- theta.names[itheta]
  if(ncol(theta.rng) > ntheta2) theta.rng <- theta.rng[,itheta,drop=FALSE]
  if(plot) {
    # set up plot
    opar <- par(no.readonly = TRUE) # save specs of current plot
    # plot size
    if(missing(layout)) {
      layout <- floor(sqrt(ntheta2))
      layout <- c(layout, ceiling(ntheta2/layout))
    }
    par(mfrow = layout, mar = c(2,2.5,2.5,0), oma = c(3, 3, .5, .5))
    on.exit(par(opar)) # restore plot parameters when exiting function
  }
  # for loop for plotting
  for(ii in 1:ntheta2) {
    ith <- itheta[ii]
    theta.seq <- seq(from = theta.rng[1,ii],
                     to = theta.rng[2,ii], len = npts)
    for(jj in 1:2) {
      # evaluate likelihood fixing all components except one
      theta.ll <- sapply(theta.seq, function(thetai) {
        theta <- theta.mle
        theta[ith] <- thetai
        loglik(theta)
      })
      if(jj == 1 && refit) {
        vth <- !is.na(theta.ll) & theta.ll > -Inf # valid values
        lth <- theta.seq < theta.mle[ith] # on the left of mle
        rth <- theta.seq > theta.mle[ith] # on the right
        # larger of the min value on each size
        lbd <- max(min(theta.ll[vth & lth]), min(theta.ll[vth & rth]))
        # rescale theta.seq to be on this range
        ibd <- c(which.min(ifelse(vth & lth, abs(theta.ll-lbd), Inf)),
                 which.min(ifelse(vth & rth, abs(theta.ll-lbd), Inf)))
        theta.seq <- seq(theta.seq[ibd[1]], theta.seq[ibd[2]], len = 200)
      } else break
    }
    # numerical check
    is.mle[ii] <- ll.max >= max(theta.ll)
    if(plot) {
      # plot loglikelihood and add estimated value
      graphics::plot(theta.seq, theta.ll, type = "l",
                     xlab = "", ylab = "")
      title(main = theta.names[ii], cex.main = 2)
      abline(v = theta.mle[ith], col = "red")
    }
  }
  if(plot) {
    # labels in margin
    mtext(side = 2, text = "Log-Likelihood",
          line = 1, outer = TRUE)
    mtext(side = 1, text = "Parameter",
          line = 1, outer = TRUE)
    return(invisible(is.mle))
  } else {
    return(is.mle)
  }
}

#' Solve method for variance matrices.
#'
#' @param V Variance matrix
#' @param x Optional vector or matrix for which to solve system of equations.  If missing calculates inverse matrix.
#' @param ldV Optionally compute log determinant as well.
#' @return Matrix solving system of equations and optionally the log-determinant.
#' @details This function is faster and more stable than \code{solve} when \code{V} is known to be positive-definite.
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V) # cholesky decomposition
  if(missing(x)) x <- diag(nrow(V))
  # solve is O(ncol(C) * ncol(x)) with triangular matrices
  # using backward subsitution
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

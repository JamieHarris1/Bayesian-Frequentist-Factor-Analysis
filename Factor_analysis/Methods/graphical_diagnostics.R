#' @title Graphical diagnostics for arrays of MCMC output.
#'
#' @description
#' This function generates graphical diagnostics for an array of MCMC output. such
#' For every parameter, a row of three plots is generated: a trace plot,
#' an ACF plot and a kernel density plot. If there is output from more than one
#' chain, the default behaviour is to overlay the plots for different chains in
#' different colours.
#'
#' @param mcmc A matrix with dimensions: iterations, parameters; or a three
#' dimensional array with dimensions: iterations, chains, parameters. The final
#' (i.e. parameter) component of the \code{dimnames} attribute of the matrix or
#' array should contain the parameter names.
#' @param rows A number indicating the number of parameters to plot per page on
#' the graphics device.
#' @param lag.max A number indicating the maximum lag for the ACF plots.
#' @param pool A logical. If \code{TRUE} the samples are pooled across chains
#' before generating the plots.
#' @param colours A vector indicating the colours to use to represent each chain.
#' Colours can be specified using any of the three kinds of R colour
#' specifications, i.e. a colour name (as listed by \code{colors()}), a
#' hexadecimal string of the form \code{"#rrggbb"} or \code{"#rrggbbaa"} or a
#' positive integer \code{i} meaning \code{palette()[i]}.
#'
#' @return \code{NULL}
#'
#' @examples
#' srs = array(rnorm(8000), c(1000, 2, 4)) # Example for illustration only!
#' dimnames(srs) = list(NULL, NULL, paste("theta[",1:4,"]",sep=""))
#' diagnostics(srs)
#'
#' @import plyr
#' @importFrom graphics lines par plot
#' @importFrom stats acf density ts var
#' @export
diagnostics = function(mcmc, rows=3, lag.max=50, pool=FALSE, colours=NULL) {
  op = par()[c("mfrow", "ask")]
  par(mfrow=c(3,3), ask=FALSE) # Avoid prompt for first plot
  if(!is.array(mcmc)) {
    stop("mcmc object must be an array.")
  }
  dims = dim(mcmc)
  if(length(dims) == 2) {
    dnames = dimnames(mcmc)
    dim(mcmc) = c(dims[1], 1, dims[2])
    dimnames(mcmc) = list(NULL, NULL, dnames[[2]])
  }
  names = dimnames(mcmc)[[3]]
  nchain = dim(mcmc)[2]
  if(is.null(colours)) {
    colours = 1:nchain
  } else {
    if(length(colours)!=nchain) {
      stop("Colour vector must have the same length as number of chains.")
    }
  }
  if(pool) {
    nm = dimnames(mcmc)
    length(nm[[2]]) = 1
    mcmc = array(mcmc, c(nchain*dim(mcmc)[1], 1, dim(mcmc)[3]), dimnames=nm)
    nchain = 1
  }
  if(nchain>1) colours_alpha = add.alpha(colours, alpha=0.4)
  else colours_alpha = colours
  count = 0
  for (name in names) {
    draws = as.matrix(mcmc[,,name])
    plot(ts(draws[,1]), main=paste("Trace plot for",name), ylab="Value",
         xlab="Iteration", ylim=c(min(draws), max(draws)), col=colours[1])
    if(nchain>1) {
      for(i in 2:nchain) {
        lines(ts(draws[,i]), col=colours[i])
      }
    }
    v = var(draws[,1])
    if (v>1e-10) {
      acfs = plyr::alply(draws, 2, acf, plot=FALSE, lag.max=lag.max)
      ymin = min(sapply(acfs, function(x) min(x$acf[,,1], na.rm = TRUE)))
      ymax = max(sapply(acfs, function(x) max(x$acf[,,1], na.rm = TRUE)))
      plot(acfs[[1]], ci=0, main=paste("ACF plot for",name), ylim=c(ymin, ymax), col=colours_alpha[1])
      if(nchain>1) {
        for(i in 2:nchain) {
          lines(0:lag.max, acfs[[i]]$acf[,,1], col=colours_alpha[i], type="h")
        }
      }
      dens = plyr::alply(draws, 2, density)
      xmin = min(sapply(dens, function(x) min(x$x)))
      xmax = max(sapply(dens, function(x) max(x$x)))
      ymin = min(sapply(dens, function(x) min(x$y)))
      ymax = max(sapply(dens, function(x) max(x$y)))
      plot(dens[[1]], main=paste("Density for",name), xlab="Value", col=colours[1],
           xlim=c(xmin, xmax), ylim=c(ymin, ymax))
      if(nchain>1) {
        for(i in 2:nchain) {
          lines(dens[[i]], col=colours[i], type="l")
        }
      }
    } else {
      plot(0,0)
      plot(0,0)
    }
    count = count + 1
    if(count==rows) par(ask=TRUE) # Now require prompt
  }
  par(op)
  NULL
}

# Internal function: add an alpha value to a colour
#' @importFrom grDevices col2rgb rgb
add.alpha = function(col, alpha=1) {
  if(missing(col)) {
    stop("Please provide a vector of colours.")
  }
  return(apply(sapply(col, col2rgb)/255, 2,
               function(x) rgb(x[1], x[2], x[3], alpha=alpha)))
}

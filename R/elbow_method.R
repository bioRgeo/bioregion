# KneeArrower: Finds Cutoff Points on Knee Curves
# Copyright 2018, 2019, 2020 Alan Tseng
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


.derivative <- function(x, y, m=0, n=50) {
  xmin <- min(x)
  xmax <- max(x)

  delta_x <- (xmax-xmin)/(n-1)
  new_x <- seq(xmin, xmax, length.out=n)
  xy <- as.data.frame(stats::approx(x, y, new_x, rule=2))
  # Do another round of smoothing
  sp <- stats::smooth.spline(xy$x, xy$y)
  xy$y <- sp$y

  new_y <- .sgolayfilt(xy$y, m=m, ts=delta_x)
  stats::approxfun(new_x, new_y)
}

.inverse <- function(f, domain) {
  function(y) {
    # Minimize |f(x)-y|
    opt <- stats::optimize(function(x) {
      abs(f(x)-y)
    }, domain)
    # If we couldn't optimize it to 0, then no solution
    # if (opt$objective < .Machine$double.eps^0.1) {
    #   opt$minimum
    # } else {
    #   NA
    # }
    opt$minimum
  }
}

.findCutoffFirstDerivative <- function(x, y, slope_ratio=0.5) {
  yf <- .derivative(x, y, 0)
  yp <- .derivative(x, y, 1)
  # Find the steepest slope either up or down
  xrange <- range(x)
  max_slope <- stats::optimize(yp, xrange, maximum=TRUE)$objective
  min_slope <- stats::optimize(yp, xrange)$objective
  steepest <- if (abs(max_slope) > abs(min_slope)) {
    max_slope
  } else {
    min_slope
  }
  # Want to find x that has the required slope
  slope <- steepest * slope_ratio
  yi <- .inverse(yp, xrange)
  knee_x <- yi(slope)
  list(x=as.numeric(knee_x), y=as.numeric(yf(knee_x)))
}

.findCutoffCurvature <- function(x, y) {
  yf <- .derivative(x, y, 0)
  yp <- .derivative(x, y, 1)
  ypp <- .derivative(x, y, 2)
  curvature <- function(x) {
    abs(ypp(x)/(1+yp(x)^2)^(3/2))
  }
  knee_x <- stats::optimize(curvature, range(x), maximum=TRUE)$maximum
  list(x=as.numeric(knee_x), y=as.numeric(yf(knee_x)))
}

.findCutoff <- function(x, y, method="first", frac.of.steepest.slope=0.5) {
  stopifnot(length(x) == length(y),
            length(x) >= 4,
            !(any(is.na(x)) || any(is.na(y)) || any(is.infinite(x)) || any(is.infinite(y))))
  if (method == "first") {
    stopifnot(0 < frac.of.steepest.slope,
              frac.of.steepest.slope <= 1)
    .findCutoffFirstDerivative(x, y, frac.of.steepest.slope)
  } else if (method == "curvature") {
    .findCutoffCurvature(x, y)
  } else {
    stop("Method must be either 'first' or 'curvature'.")
  }
}


# Functions below come from package signal
# License GPL-V2
# Original Octave version by Paul Kienzle pkienzle@users.sf.net. Modified by Pascal Dupuis. Conversion to R by Tom Short.
.sgolayfilt <- function (x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1)
{
  len = length(x)
  if (class(p) == "sgolayFilter" || (!is.null(dim(p)) &&
                                     dim(p) > 1)) {
    F = p
    n = nrow(F)
  }
  else F = .sgolay(p, n, m, ts)
  k = floor(n/2)
  z = .filter(F[k + 1, n:1], 1, x)
  c(F[1:k, ] %*% x[1:n], z[n:len], F[(k + 2):n, ] %*% x[(len -
                                                           n + 1):len])
}

.sgolay <- function (p, n, m = 0, ts = 1)
{
  if (n%%2 != 1)
    stop("sgolay needs an odd filter length n")
  if (p >= n)
    stop("sgolay needs filter length n larger than polynomial order p")
  Fm <- matrix(0, n, n)
  k <- floor(n/2)
  for (row in 1:(k + 1)) {
    Ce <- (((1:n) - row) %*% matrix(1, 1, p + 1))^(matrix(1,
                                                          n) %*% (0:p))
    A <- MASS::ginv(Ce, tol = .Machine$double.eps)
    Fm[row, ] <- A[1 + m, ]
  }
  Fm[(k + 2):n, ] <- (-1)^m * Fm[k:1, n:1]
  if (m > 0)
    Fm <- Fm * prod(1:m)/(ts^m)
  class(Fm) <- "sgolayFilter"
  Fm
}


.filter <- function(filt, a, x, init, init.x, init.y, ...) {
  if(missing(init.x))
    init.x <- c(rep(0, length(filt) - 1))
  if(length(init.x) != length(filt) - 1)
    stop("length of init.x should match filter length-1 = ", length(filt)-1)
  if(missing(init) && !missing(init.y))
    init <- rev(init.y)
  if(all(is.na(x)))
    return(x)
  if (length(filt)) {
    x1 <- stats::filter(c(init.x, x), filt / a[1], sides = 1)
    if(all(is.na(x1)))
      return(x)
    x <- stats::na.omit(x1, filt / a[1] , sides = 1)
  }
  if (length(a) >= 2)
    x <- stats::filter(x, -a[-1] / a[1], method = "recursive", init = init)
  x
}

#' Wrap an angle into [0, 2*pi)
#'
#' @title Adjust angle to [0, 2*pi)
#'
#' @description
#' Normalises a single numeric angle in radians to the half-open interval
#' \[0, 2*pi) by adding or subtracting the appropriate multiple of 2*pi.
#' Used internally after \code{atan2} to convert phase estimates to the
#' standard Zeitgeber-time scale.
#'
#' @param x Numeric scalar; angle in radians.
#'
#' @return Numeric scalar in \[0, 2*pi).
#'
#' @export
#'
#' @examples
#' adjust.to.2pi(-pi/2)  # returns 3*pi/2
#' adjust.to.2pi(5*pi)   # returns pi
adjust.to.2pi = function(x){
  ((x %% (2*pi)) + 2*pi) %% (2*pi)
}

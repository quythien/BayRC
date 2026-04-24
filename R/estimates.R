#' Extract posterior estimates and credible intervals for all parameters
#'
#' @title Get all posterior estimates from BayRC MCMC output
#'
#' @description
#' Computes posterior estimates and highest-density credible intervals for
#' amplitude (A), acrophase (phi), MESOR (M), residual variance (sigma),
#' and optionally peak time (t_p) from the MCMC output list.  Circular
#' statistics (\code{circular_median} and \code{circular_HDI}) are used
#' for phi/t_p; linear posterior means and HDI are used for the rest.
#'
#' @param res Named list; MCMC output from \code{CB_MCMC_single_rj_slice}
#'   (or \code{CBt_MCMC_single}).  Must contain matrices \code{rho},
#'   \code{A}, \code{phi}, \code{M}, \code{sigma} with
#'   \code{attr(rho, "symbols")} and \code{attr(rho, "RHYindex")} set by
#'   \code{match_symbols}.
#' @param burn Integer; number of initial MCMC iterations to discard as
#'   additional burn-in when computing estimates (default 100).
#' @param CI Logical; if \code{TRUE} (default), return lower and upper
#'   credible-interval bounds in addition to the point estimate.
#' @param credMass Numeric in (0, 1); nominal coverage of the HDI
#'   (default 0.95).
#' @param P Numeric; period in hours used for circular statistics
#'   (default 24).
#' @param rhythmic Logical; if \code{TRUE}, return only rows for genes
#'   classified as rhythmic (RHYindex == 1; default \code{FALSE}).
#'
#' @return A list of data.frames, one per parameter, each with columns
#'   \code{<param>.Est}, \code{<param>.Lower}, \code{<param>.Upper} and
#'   \code{RHYindex}.  Row names are gene symbols when
#'   \code{attr(rho, "symbols")} is set.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim  <- CBt_sim_data()
#' dat  <- list(data  = as.data.frame(sim$data[[1]]$dat),
#'              time  = sim$data[[1]]$x$time,
#'              gname = paste0("G", seq_len(nrow(sim$data[[1]]$dat))))
#' init <- CB_init_single(dat)
#' res  <- CB_MCMC_single_rj_slice(dat, init,
#'           iteration = 300, thin = 10, n.burn = 100)
#' est  <- CB_getAllEst(res, burn = 10)
#' }
CB_getAllEst <- function(res, burn=100, CI=TRUE, credMass=0.95, P=24, rhythmic=FALSE) {
  params <- c("A", "phi", "M", "sigma")
  if ("t_p" %in% names(res)) params <- c(params, "t_p")
  
  # Extract symbols and RHYindex once
  symbols <- attr(res$rho, "symbols")
  RHYindex <- attr(res$rho, "RHYindex")
  
  res.all <- lapply(params, function(a.param) {
    a.MCMC <- res[[a.param]]
    
    # **Compute estimates for all rows first**
    a.tab <- t(sapply(1:nrow(a.MCMC), function(i) {
      get_est_CI(a.MCMC[i, ], res$rho[i, ], burn, a.param, credMass, P)
    }))
    
    # Convert to DataFrame and set column names
    a.tab <- as.data.frame(a.tab)
    colnames(a.tab) <- paste0(a.param, c(".Est", ".Lower", ".Upper"))
    
    # Assign row names
    if (!is.null(symbols)) rownames(a.tab) <- symbols
    
    # Add RHYindex column
    a.tab$RHYindex <- RHYindex
    
    # **Only filter at the end** (reducing memory allocation)
    if (rhythmic) a.tab <- a.tab[RHYindex == 1, , drop=FALSE]
    
    if (!CI) a.tab <- a.tab[, 1, drop=FALSE]
    
    return(a.tab)
  })
  
  return(res.all)
}

#' Compute posterior estimate and HDI for one parameter row
#'
#' @description
#' Dispatcher that routes to the circular estimator for phi/t_p or the
#' linear estimator for A/M/sigma, applying burn-in trimming and
#' rhythmicity masking for amplitude.
#'
#' @param MCMC Numeric vector; posterior samples for one gene.
#' @param rho Numeric vector; posterior rhythmicity indicators for the same gene.
#' @param burn Integer; additional burn-in samples to discard.
#' @param a.param Character; parameter name (\code{"A"}, \code{"phi"},
#'   \code{"M"}, \code{"sigma"}, or \code{"t_p"}).
#' @param credMass Numeric; HDI coverage.
#' @param P Numeric; period for circular wrapping.
#'
#' @return Numeric vector of length 3: estimate, lower CI, upper CI.
#'
#' @keywords internal
get_est_CI = function(MCMC, rho, burn, a.param, credMass, P){
  # t = False for peak time 
  if(a.param=="phi"){
    out = get_t_phi_CI_est(MCMC, P, credMass, burn, rho, t=FALSE)
  }else if (a.param=="t_p"){
    out = get_t_phi_CI_est(MCMC, P, credMass, burn, rho=1, t=TRUE)
  }else{
    MCMCsub = MCMC[-seq_len(burn)]
    if(a.param=="A"){
      rhosub = rho[-seq_len(burn)]
      MCMCsub = MCMCsub[rhosub==1]
    }
    out = get_other_CI_est(MCMCsub, credMass)
  }
  return(out)
}

#' Compute mean and HDI for a linear parameter
#'
#' @description Returns the posterior mean and highest-density interval for
#' a non-circular parameter using \code{HDInterval::hdi}.
#'
#' @param MCMC Numeric vector of posterior samples.
#' @param credMass Numeric; HDI coverage probability.
#'
#' @return Numeric vector of length 3: mean, lower HDI bound, upper HDI bound.
#'
#' @keywords internal
get_other_CI_est = function(MCMC, credMass){
  est = mean(MCMC, na.rm=TRUE)
  CI = HDInterval::hdi(MCMC, credMass)
  c(est, CI)
}

#' Compute the circular median of phase samples
#'
#' @description
#' Converts samples from hours to radians, computes the circular median
#' via \code{circular::median.circular}, and converts back to hours.
#'
#' @param samples Numeric vector; phase samples in hours.
#' @param P Numeric; period in hours (default 24).
#'
#' @return Numeric scalar; circular median in hours in \[0, P).
#'
#' @export
circular_median <- function(samples, P = 24) {
  samples_rad <- circular((samples / P) * (2 * pi), units = "radians")
  median_rad <- median.circular(samples_rad, type = "median")
  median_time <- as.numeric(median_rad) / (2 * pi) * P
  return(median_time %% P)
}

#' Normalize an angle to a specified interval
#'
#' @description
#' Wraps a numeric angle \code{x} into the half-open interval
#' \[a, a + P) using modular arithmetic.  Useful for mapping acrophase
#' samples to a canonical range such as \[0, 24) or \[-12, 12).
#'
#' @param x Numeric; angle(s) to normalise (hours).
#' @param a Numeric; left endpoint of the target interval (default 0).
#' @param P Numeric; period / interval width (default 24).
#'
#' @return Numeric; normalised angle(s) in \[a, a + P).
#'
#' @export
normalize_angle <- function(x, a = 0, P = 24) {
  return(((x - a) %% P) + a)
}

#' Recenter phase samples around their circular median
#'
#' @description
#' Shifts samples so that their circular median is at zero, then wraps
#' into \[a, a + P).  Used internally before the sliding-window HDI.
#'
#' @param samples Numeric vector; phase samples in hours.
#' @param P Numeric; period (default 24).
#' @param a Numeric; left endpoint of normalisation interval (default 0).
#'
#' @return List with \code{shifted} (re-centred samples) and \code{med}
#'   (the circular median used for re-centring).
#'
#' @keywords internal
recenter_samples <- function(samples, P = 24, a = 0) {
  samples_rad <- (samples / P) * (2 * pi)
  median_rad <- median.circular(circular(samples_rad), type = "median")
  med <- (median_rad / (2 * pi)) * P
  shifted <- samples - med
  shifted <- normalize_angle(shifted, a = a, P = P)
  return(list(shifted = shifted, med = med))
}

#' Compute circular estimate and HDI for phase or peak-time
#'
#' @description
#' Returns the circular median and \code{circular_HDI} bounds for a single
#' gene's phase (phi) or peak-time (t_p) MCMC chain, after discarding
#' burn-in and optionally filtering by rhythmicity.
#'
#' @param row_MCMC Numeric vector; posterior samples for one gene.
#' @param P Numeric; period (default 24).
#' @param credMass Numeric; HDI coverage (default 0.95).
#' @param burn Integer; burn-in samples to discard (default 1).
#' @param rho Numeric vector; rhythmicity samples (used to gate A samples
#'   when \code{t = FALSE}).  Set to 1 for t_p.
#' @param t Logical; if \code{TRUE} compute for peak time; if \code{FALSE}
#'   compute for acrophase phi.
#' @param a Numeric; left endpoint of the normalisation interval (default 0).
#'
#' @return Named numeric vector with elements \code{phi.Est},
#'   \code{phi.Lower}, \code{phi.Upper}.
#'
#' @keywords internal
get_t_phi_CI_est <- function(row_MCMC, P = 24, credMass = 0.95, burn = 1, rho = 1, t = TRUE, a = 0) {
  
  if (length(row_MCMC) <= burn) {
    stop("row_MCMC has insufficient data after burn-in. Check input.")
  }
  row_MCMC <- row_MCMC[-seq_len(burn)]

  phi.Est <- circular_median(row_MCMC, P = P)
  hdi_bounds <- circular_HDI(row_MCMC, credMass = credMass, P = P, a = a)
  phi.Lower <- hdi_bounds$lower
  phi.Upper <- hdi_bounds$upper

  if (phi.Upper < phi.Lower) {
    phi.Upper <- normalize_angle(phi.Upper + P, a = a, P = P)
  }

  return(c(phi.Est = phi.Est, phi.Lower = phi.Lower, phi.Upper = phi.Upper))
}

# Old
# circular_HDI <- function(samples, credMass = 0.95, P = 24, a = 0) {
#   # Recenter the samples and extract the median.
#   rec <- recenter_samples(samples, P = P, a = a)
#   shifted <- rec$shifted
#   med <- rec$med
#   
#   # Sort the shifted samples.
#   sorted_samples <- sort(as.numeric(shifted))
#   n <- length(sorted_samples)
#   interval_size <- floor(credMass * n)
#   
#   # Create an extended array to capture wrap-around.
#   extended_samples <- c(sorted_samples, sorted_samples + P)
#   
#   min_width <- Inf
#   best_lower <- NA
#   best_upper <- NA
#   
#   # Slide a window of size interval_size over the extended array.
#   for (i in 1:n) {
#     j <- i + interval_size
#     if (j <= length(extended_samples)) {
#       width <- extended_samples[j] - extended_samples[i]
#       if (width < min_width) {
#         min_width <- width
#         best_lower <- extended_samples[i]
#         best_upper <- extended_samples[j]
#       }
#     }
#   }
#   
#   # Shift the interval back by adding the median.
#   CI_lower <- best_lower + med
#   CI_upper <- best_upper + med
#   
#   # Normalize back to [a, a+P)
#   CI_lower <- normalize_angle(CI_lower, a = a, P = P)
#   CI_upper <- normalize_angle(CI_upper, a = a, P = P)
#   
#   return(list(lower = CI_lower, upper = CI_upper))
#   
# }

# New

#' Compute the shortest circular highest-density interval
#'
#' @description
#' Finds the shortest interval on a circle of circumference \code{P} that
#' contains at least \code{credMass} of the posterior samples.  Uses a
#' sliding-window algorithm on a doubled sorted array to handle wrap-around.
#'
#' @param samples Numeric vector; phase samples in hours.
#' @param credMass Numeric; desired coverage probability (default 0.95).
#' @param P Numeric; period / circle circumference in hours (default 24).
#' @param a Numeric; left endpoint of the normalisation interval (default 0).
#'
#' @return List with elements \code{lower} and \code{upper} (both in
#'   \[a, a + P)) giving the bounds of the shortest HDI.
#'
#' @export
circular_HDI <- function(samples, credMass = 0.95, P = 24, a = 0) {
  # Step 1: wrap samples into [a, a+P)
  samples <- normalize_angle(samples, a = a, P = P)
  sorted  <- sort(samples)
  n <- length(sorted)
  interval_size <- floor(credMass * n)
  
  # Step 2: duplicate the circle to handle wrap-around
  extended <- c(sorted, sorted + P)
  
  # Step 3: sliding window to find shortest interval
  min_width <- Inf
  best_lower <- NA
  best_upper <- NA
  
  for (i in 1:n) {
    j <- i + interval_size
    if (j <= length(extended)) {
      width <- extended[j] - extended[i]
      if (width < min_width) {
        min_width <- width
        best_lower <- extended[i]
        best_upper <- extended[j]
      }
    }
  }
  
  # Step 4: normalize the result
  CI_lower <- normalize_angle(best_lower, a = a, P = P)
  CI_upper <- normalize_angle(best_upper, a = a, P = P)
  
  return(list(lower = CI_lower, upper = CI_upper))
}



#' Test whether a point lies within a circular interval
#'
#' @description
#' Returns \code{TRUE} if \code{point} lies within the circular interval
#' \[lower, upper\] on a circle of circumference \code{P}, handling
#' wrap-around (i.e. lower > upper means the interval crosses zero).
#'
#' @param lower Numeric; lower bound of the interval.
#' @param upper Numeric; upper bound of the interval.
#' @param point Numeric; test point (default 0).
#' @param P Numeric; period (default 24).
#' @param a Numeric; left endpoint of normalisation interval (default 0).
#'
#' @return Logical scalar.
#'
#' @keywords internal
in_circular_interval <- function(lower, upper, point = 0, P = 24, a = 0) {
  lower <- normalize_angle(lower, a = a, P = P)
  upper <- normalize_angle(upper, a = a, P = P)
  point <- normalize_angle(point, a = a, P = P)

  if (lower <= upper) {
    return(point >= lower && point <= upper)
  } else {
    return(point >= lower || point <= upper)
  }
}
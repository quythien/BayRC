#' @param MCMC 
#' @param burn 
#' @param credMass 
#' @param P 
#' @param phi_adj 
#'
#' @return
#' @export
#'
#' @examples
#' 



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

get_other_CI_est = function(MCMC, credMass){
  est = mean(MCMC, na.rm=TRUE)
  CI = HDInterval::hdi(MCMC, credMass)
  c(est, CI)
}

circular_median <- function(samples, P = 24) {
  samples_rad <- circular((samples / P) * (2 * pi), units="radians")
  median_rad <- median.circular(samples_rad, type = "median")
  median_time <- (median_rad / (2 * pi)) * P
  return(median_time %% P)
}

# Helper function to normalize an angle to a specified interval.
normalize_angle <- function(x, a = 0, P = 24) {
  return(((x - a) %% P) + a)
}

recenter_samples <- function(samples, P = 24, a = 0) {
  samples_rad <- (samples / P) * (2 * pi)
  median_rad <- median.circular(circular(samples_rad), type = "median")
  med <- (median_rad / (2 * pi)) * P
  shifted <- samples - med
  shifted <- normalize_angle(shifted, a = a, P = P)
  return(list(shifted = shifted, med = med))
}

get_t_phi_CI_est <- function(row_MCMC, P = 24, credMass = 0.95, burn = 1, rho = 1, t = TRUE, a = 0) {
  
  circular_median <- function(samples, P = 24) {
    samples_rad <- circular((samples / P) * (2 * pi), units="radians")
    median_rad <- median.circular(samples_rad, type = "median")
    median_time <- (median_rad / (2 * pi)) * P
    return(median_time %% P)
  }
  
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
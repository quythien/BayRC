#' Initialize MCMC chain for a dataset with Zeitgeber-time uncertainty
#'
#' @title Initialize time-error BayRC MCMC chain
#'
#' @description
#' Produces starting values for the CBt MCMC sampler, which models
#' Zeitgeber-time measurement error (t_p).  Parameter initialisation
#' follows the same OLS cosinor approach as \code{CB_init_single} but
#' also initialises the time-error latent variable \code{t_p} to zero.
#'
#' @param Data.list A named list with three elements: \code{data} (G x N
#'   data.frame of expression values), \code{time} (length-N numeric
#'   Zeitgeber time vector in hours), and \code{gname} (length-G character
#'   gene identifiers).
#' @param P Numeric; circadian period in hours (default 24).
#' @param FitCosinor Logical; if \code{TRUE} (default) use OLS cosinor
#'   estimates; if \code{FALSE} draw from priors.
#' @param mu_M Numeric; prior mean for MESOR (default 0).
#' @param sigma_M Numeric; prior standard deviation for MESOR (default 10).
#' @param mu_A Numeric; prior mean for amplitude (default 1).
#' @param sigma_A Numeric; prior standard deviation for amplitude
#'   (default 10).
#' @param seed Integer; random seed (default 15213).
#'
#' @return A named list with elements \code{rho}, \code{M}, \code{A},
#'   \code{phi}, \code{sigma} (same as \code{CB_init_single}), plus
#'   \code{t_p} (length-N vector of initial time-error values, set to 0).
#'   Suitable for passing to \code{CBt_MCMC_single} as \code{Init.value}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim  <- CBt_sim_data()
#' dat  <- list(data  = as.data.frame(sim$data[[1]]$dat),
#'              time  = sim$data[[1]]$x$time,
#'              gname = paste0("G", seq_len(nrow(sim$data[[1]]$dat))))
#' init <- CBt_init_single(dat)
#' }
CBt_init_single = function(Data.list = list(data = as.data.frame(a.sim.dat$data[[1]]$dat),
                                           time = a.sim.dat$data[[1]]$x$time,
                                           gname = paste0("G", seq_len(nrow(a.sim.dat$data[[1]]$dat)))), 
                    P = 24, FitCosinor = TRUE, 
                    mu_M = 0, sigma_M = 10, mu_A = 1, sigma_A = 10, #theta = pi, kappa = 0.01,
#                   VM_theta0 = 0, VM_R0 = 1,VM_c = 1, VM_Bessel_k = 10, 
#                   nClass.init = 10,
                    seed = 15213){
  # P = 24; mu_M = 0; sigma_M = 10; mu_A = 0; sigma_A = 10; theta = 0; kappa = 0.001;
  set.seed(seed)
  #later: data input check point
  omega = 2*pi/P
  Y = Data.list[[1]]
  t = Data.list[[2]]
  G = nrow(Y)
  N = ncol(Y)
  t.c = cos(omega*t)
  t.s = sin(omega*t)
  
  if(FitCosinor){
    a.x = cbind(1, t.c, t.s)
    xx.inv = solve(t(a.x)%*%a.x)
    a.beta = solve(t(a.x)%*%a.x)%*%t(a.x)%*%t(Y)
    out.M = a.beta[1, ]
    out.Acosphi = a.beta[2, ]
    out.Asinphi = a.beta[3, ]
    out.A = sqrt(a.beta[2, ]^2+a.beta[3, ]^2)
    out.cosphi = out.Acosphi/out.A
    out.sinphi = out.Asinphi/out.A
    out.tanphi = out.sinphi/out.cosphi
    out.phi0 = atan(out.tanphi)
    out.phi = ifelse(out.cosphi<0, out.phi0+pi, out.phi0+2*pi*out.tanphi<0) #adjust the atan output to [0, 2*pi)
    est.Y = t(a.x%*%a.beta)
    out.sigma = apply((Y - est.Y)^2, 1, mean)
    a.F = apply(a.beta[2:3, ], 2, function(a){
      a = matrix(a, ncol = 1)
      t(a)%*%solve(xx.inv[2:3, 2:3])%*%a
      })/2/out.sigma
    out.rho = pf(a.F, 2, N-3, lower.tail = FALSE)<0.05

  }else{
    out.M = rnorm(G, mu_M, sigma_M) 
    out.A = truncnorm::rtruncnorm(G, 0, 3, mu_A, sigma_A)
    # phi = matrix(Rfast::rvonmises(G*J, m = theta, k = kappa), byrow = FALSE, nrow = G)/(2*pi)*P #uniform
    out.phi = runif(G,0, P) #uniform
    out.sigma = rep(1, G) ##1. sigma is indepedent for each gene and each group
    out.rho = rbinom(G, 1, 0.5)
  }
  t_p = rep(0, N)

  init=list(rho = out.rho,
            M = out.M,
            A = out.A,
            phi = out.phi,
            sigma = out.sigma,
            t_p = t_p)
  
  print("MCMC chain initialized")
  return(init)
}
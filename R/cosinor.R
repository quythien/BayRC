Cosinor_fit = function(x, period = 24, amp.cutoff = 0,  p.adjust.method = "BH", parallel.ncores  = 1){

  data = as.matrix(x$data)
  time = x$time
  gname = x$gname
  
  G = nrow(data)
  res = lapply(1:G, function(g){
    one_cosinor_OLS(time, y=data[g, ], period)
  })
  res = do.call(rbind.data.frame, res)
  res = cbind.data.frame(gname = gname, res)
  res$qvalue = stats::p.adjust(res$pvalue, p.adjust.method)
  
  x$rhythm = res
  x$P = period
  return(x)
}

one_cosinor_OLS = function(tod = time, y = y, period = 24){
  
  n = length(tod)
  omega = 2*pi/period
  x1 = cos(omega*tod)
  x2 = sin(omega*tod)
  # mat.X = matrix(c(rep(1, n), x1, x2), ncol = 3, byrow = FALSE)
  # mat.XX = t(mat.X)%*%mat.X#mat.XX = mat.S
  mat.S = matrix(c(n, sum(x1), sum(x2),
                   sum(x1), sum(x1^2), sum(x1*x2),
                   sum(x2), sum(x1*x2), sum(x2^2)),
                 nrow = 3, byrow = TRUE)
  vec.d = c(sum(y), sum(y*x1), sum(y*x2))
  
  mat.S.inv = solve(mat.S)
  est = mat.S.inv%*%vec.d
  m.hat = est[1]
  beta1.hat = est[2]
  beta2.hat = est[3]
  A.hat = sqrt(beta1.hat^2 + beta2.hat^2)
  phase.hat = adjust.to.2pi(atan2(beta2.hat, beta1.hat))/omega

  #inference
  TSS = sum((y-mean(y))^2)
  yhat = m.hat + beta1.hat*x1+beta2.hat*x2
  RSS = sum((y-yhat)^2)
  MSS = TSS-RSS
  Fstat = (MSS/2)/(RSS/(n-3))
  pval = stats::pf(Fstat, 2, n-3, lower.tail = FALSE)
  sigma2.hat = RSS/(n-3)
  sigma.hat = sqrt(sigma2.hat)
  
  #output
  out = list(M = m.hat,
             A = A.hat, 
             phase = phase.hat, 
             peak = phase.hat, 
             R2 = MSS/TSS, 
             sigma = sigma.hat,
             pvalue = pval)
  return(out)
}

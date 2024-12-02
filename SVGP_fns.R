library(ggplot2)
library(tidyverse)
library(gdata)
library(gridExtra)
library(MASS)
library(latex2exp)
library(gdata)
library(hrbrthemes)
library(patchwork)
library(xtable)
library(data.table)
# investigating the difference between the variational posterior mean
# and the full posterior mean
#### HELPER FUNCTIONS ####
expand.grid.extra = function(params_1, params_2){
  #expand.grid, but adds each case of params_2 to each case of params_1
  final_params = matrix(nrow = dim(params_1)[1]*dim(params_2)[1],
                        ncol = dim(params_1)[2] + dim(params_2)[2])
  for(i in 1:dim(params_1)[1]){
    for(j in 1:dim(params_2)[1]){
      param_ix = (i-1)*dim(params_2)[1] + j
      # print(c(as.matrix(params_1[i,]), as.matrix(params_2[j,])))
      final_params[param_ix,] = c(as.matrix(params_1[i,]), as.matrix(params_2[j,]))
    }
  }
  colnames(final_params) = c(colnames(params_1), colnames(params_2))
  return(final_params)
}
#### KERNELS ####
bm_kern = function(x, y, cn = 1){
  #rescaled Brownian motion kernel
  return(cn*min(x, y))
}
mat_kern_3_2 = function(x, y, cn = 1){
  r = sqrt(sum((x - y)^2))
  return( (1+sqrt(3)*r/cn)*exp(-sqrt(3)*r/cn) )
}
mat_kern_5_2 = function(x, y, cn = 1){
  r = sqrt(sum((x - y)^2))
  return( (1+sqrt(5)*r/cn+ 5*r^2/(3*cn^2))*exp(-sqrt(5)*r/cn) )
}
mat_kern_1_2 = function(x, y, cn = 1){
  r = sqrt(sum((x - y)^2))
  return(exp(-r/cn))
}
se_kern = function(x, y, cn = 1){
  # provide cn as the lengthscale
  return(exp(-(x-y)^2 / cn^2))
}

se_kern_multi_d = function(x, y, cn = 1){
  # provide cn as the lengthscale
  r = sqrt(sum((x - y)^2))
  return(exp(-r^2 / cn^2))
}

lin_kern = function(x, y){
  return( 1+ x*y ) 
}

#### BASIS FUNCTIONS ####
basis_functions_x = function(K, xs){
  # K is the (even) number of basis functions to compute 
  args = seq(1, K/2, by=1)*2*pi
  args = as.matrix(xs)%*%t(as.matrix(args))
  # L is the size of the range discretisation
  cosines = sqrt(2)*cos(args)
  sines = sqrt(2)*sin(args)
  ones = as.matrix(rep(1,length(xs)))
  X = cbind(ones, t(interleave(t(cosines), t(sines))))
  return(X)
}
compute_fn_from_coefs_x = function(coefs, xs){
  K = length(coefs) - 1
  basis_fn_vals = basis_functions_x(K,xs)
  fn_vals = basis_fn_vals%*%coefs
  return(fn_vals)
}
f0_xs = function(xs, alpha){
  if(alpha <= 1){
    return(abs(xs - 1/2)^alpha)
  }else if(alpha == 2){
    return( sign(xs - 1/2) * abs(xs - 1/2)^alpha_2   )
  }
}

# f0_xs = function(xs, alpha){
#   if(alpha == 1){
#     return(sign(xs - 1/2)*(abs(xs - 1/2)^2) )
#   }else if(alpha == 2){
#     return( abs(xs - 1/2)^3)
#   }
# }

f0_xs_multi_d = function(Xs, alpha){
  if(!is.matrix(Xs)){
    Xs = matrix(Xs, nrow = 1)
  }
  apply(Xs, 1, function(x) sqrt(sum(x^2))^alpha )
}

#### VARIATIONAL POSTERIOR FUNCTIONS ####

vp_mean_variance = function(x_, m_, svd_, k_, xns_, y_, cn_ = 1, sigma=1){
  eta_ks = 1/(svd_$d[1:m_] + sigma^2)
  V_m = svd_$u[,1:m_]
  A = V_m %*% diag(eta_ks) %*% t(V_m)
  k_n_x = sapply(xns_, function(y) k_(y, x_, cn = cn_))
  mean_x = t(k_n_x) %*% A %*% y_
  var_x = k_(x_, x_, cn = cn_) - t(k_n_x) %*% A %*% k_n_x
  return(c(mean_x, var_x))
}

vp_mean_variance_multiple_m = function(x_, m_s, svd_, k_, xns_, y_, cn_ = 1, sigma=1){
  ret_mat = matrix(ncol = 2, nrow = length(m_s))
  k_n_x = sapply(xns_, function(y) k_(y, x_, cn = cn_))
  for(ix in 1:length(m_s)){
    m_ = m_s[ix]
    eta_ks = 1/(svd_$d[1:m_] + sigma^2)
    V_m = svd_$u[,1:m_]
    A = V_m %*% diag(eta_ks) %*% t(V_m)
    mean_x = t(k_n_x) %*% A %*% y_
    var_x = k_(x_, x_, cn = cn_) - t(k_n_x) %*% A %*% k_n_x
    ret_mat[ix,1] = mean_x
    ret_mat[ix,2] = var_x
  }
  return(ret_mat)  
}

vp_mean_variance_multiple_m_multi_d = function(x_, m_s, svd_, k_, xns_, y_, cn_ = 1, sigma = 1){
  ret_mat = matrix(ncol = 2, nrow = length(m_s))
  k_n_x = apply(xns_, 1, function(x) k_(x, x_, cn = cn_))
  for(ix in 1:length(m_s)){
    m_ = m_s[ix]
    eta_ks = 1/(svd_$d[1:m_] + sigma^2)
    V_m = svd_$u[,1:m_]
    A = V_m %*% diag(eta_ks) %*% t(V_m)
    mean_x = t(k_n_x) %*% A %*% y_
    var_x = k_(x_, x_, cn = cn_) - t(k_n_x) %*% A %*% k_n_x
    ret_mat[ix,1] = mean_x
    ret_mat[ix,2] = var_x
  }
  return(ret_mat)  
}

sample_variational_CS = function(gamma_, m_, x0_, svd_, alpha_, cn_ = 1, rand_design=FALSE){
  n = length(svd_$d)
  if(rand_design){
    xns = runif(n, 0, 1) 
  }else{
    xns = c(1:n)/(n+0.5)  
  }
  K = 401
  #coeffs of f0
  f0_ks = (1:K)^(-1/2 - alpha)
  f0 = compute_fn_from_coefs_x(f0_ks, xs = xns)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  
  sigma_est = maximise_sigma_known_svd(y, svd_)
  
  vp_fit = vp_mean_variance(x_ = x0_, m_=m_, svd_=svd_, k_=kern, xns_=xns,
                            y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[1]
  vp_var = vp_fit[2]
  
  z_gamma = qnorm((1+gamma)/2)
  vp_CS = c(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
  # f0_x0 = compute_fn_from_coefs_x(f0_ks, xs = x0_)
  # 
  # # see if f0_x0 is in vp_CS
  # hit = vp_CS[1] <= f0_x0 & f0_x0 <= vp_CS[2]
  # return()
}

sample_variational_CS_multiple_m = function(gamma_, m_s, x0_, svd_, alpha_, cn_ = 1, kern = bm_kern,
                                            xns, noise = 'gaussian'){
  n = dim(svd_$u)[1]
  # if(rand_design){
  #   xns = runif(n, 0, 1) 
  # }else{
  #   xns = c(1:n)/(n+0.5)  
  # }
  # K = 401
  #coeffs of f0
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0 = compute_fn_from_coefs_x(f0_ks, xs = xns)
  f0 = f0_xs(xs = xns, alpha = alpha_)
  if(noise == 'gaussian'){
    eps = rnorm(n)
  }else if(noise == 'laplace'){
    eps = rlaplace(n)
  }
  y = f0 + eps
  sigma_est = maximise_sigma_known_svd(y, svd_)
  cat('Estimated sigma: ', sigma_est, '\n')
  vp_fit = vp_mean_variance_multiple_m(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                       y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  vp_CS = cbind(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
}

compute_NLPD = function(gamma_, m_s, x0_, svd_, alpha_, cn_ = 1, kern = bm_kern,
                        xns, noise = 'gaussian'){
  n = dim(svd_$u)[1]
  f0 = f0_xs(xs = xns, alpha = alpha_)
  if(noise == 'gaussian'){
    eps = rnorm(n)
  }else if(noise == 'laplace'){
    eps = rlaplace(n)
  }
  y = f0 + eps
  sigma_est = maximise_sigma_known_svd(y, svd_)
  vp_fit = vp_mean_variance_multiple_m(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                       y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  NLPDs = -log(dnorm(0, vp_mean, sd = sqrt(vp_var)))
  return(NLPDs)
}


sample_variational_CS_multiple_m_random_design = function(gamma_, m_s, x0_, n, alpha_, cn_ = 1, kern = bm_kern, noise = 'gaussian'){
  xns = runif(n)
  f0 = f0_xs(xs = xns, alpha = alpha_)
  if(noise == 'gaussian'){
    eps = rnorm(n)
  }else if(noise == 'laplace'){
    eps = rlaplace(n)
  }
  y = f0 + eps
  K_nn = sapply(xns, function(x) sapply(xns, function(y) kern(x, y, cn = cn_)))
  svd_ = svd(K_nn)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  vp_fit = vp_mean_variance_multiple_m(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                       y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  vp_CS = cbind(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
}

compute_NLPD_random_design = function(gamma_, m_s, x0_, n, alpha_, cn_ = 1, kern = bm_kern, noise = 'gaussian'){
  xns = runif(n)
  f0 = f0_xs(xs = xns, alpha = alpha_)
  if(noise == 'gaussian'){
    eps = rnorm(n)
  }else if(noise == 'laplace'){
    eps = rlaplace(n)
  }
  y = f0 + eps
  K_nn = sapply(xns, function(x) sapply(xns, function(y) kern(x, y, cn = cn_)))
  svd_ = svd(K_nn)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  vp_fit = vp_mean_variance_multiple_m(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                       y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  NLPDs = -log(dnorm(0, vp_mean, sd = sqrt(vp_var)))
  return(NLPDs)
}

sample_variational_CS_multiple_m_multi_d = function(gamma_, m_s, x0_, svd_, alpha_, cn_ = 1, kern = bm_kern,
                                                    xns){
  n = dim(svd_$u)[1]
  # K = 401
  # #coeffs of f0
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0 = apply(xns, 1, function(x) sum(sapply(x, function(xi) compute_fn_from_coefs_x(f0_ks, xs = xi))))
  f0 = f0_xs_multi_d(Xs = xns, alpha = alpha_)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  # sigma_est = maximise_sigma_known_svd(y, svd_)
  sigma_est = 1
  
  vp_fit = vp_mean_variance_multiple_m_multi_d(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                               y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  vp_CS = cbind(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
}
sample_variational_CS_multiple_m_multi_d_adaptive = function(gamma_, m_s, x0_, alpha_, kern = mat_kern_1_2,
                                                             xns){
  n = dim(xns)[1]
  # K = 401
  # #coeffs of f0
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0 = apply(xns, 1, function(x) sum(sapply(x, function(xi) compute_fn_from_coefs_x(f0_ks, xs = xi))))
  f0 = f0_xs_multi_d(Xs = xns, alpha = alpha_)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  cn_est = maximise_cn(y, xns, kern = kern)
  # cn_est = n^{-1/(1+2*alpha_)}
  
  K_nn = apply(xns, 1, function(x) apply(xns, 1, function(y) kern(x, y, cn)))
  svd_ = svd(K_nn)
  
  vp_fit = vp_mean_variance_multiple_m_multi_d(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                               y_=y, cn_ = cn_est, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  vp_CS = cbind(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
}

compute_NLPD_multi_d = function(gamma_, m_s, x0_, svd_, alpha_, cn_ = 1, kern = bm_kern,
                                xns){
  n = dim(svd_$u)[1]
  # K = 401
  # #coeffs of f0
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0 = apply(xns, 1, function(x) sum(sapply(x, function(xi) compute_fn_from_coefs_x(f0_ks, xs = xi))))
  f0 = f0_xs_multi_d(Xs = xns, alpha = alpha_)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  
  vp_fit = vp_mean_variance_multiple_m_multi_d(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                               y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  NLPDs = -log(dnorm(0, vp_mean, sd = sqrt(vp_var)))
  return(NLPDs)
}

sample_correlated_inputs = function(n, p, rho){
  print(rho)
  Sigma = matrix(rho, nrow = p, ncol = p)
  diag(Sigma) = rep(1,p)
  U = chol(Sigma)
  normal_sample = matrix(rnorm(n*p), nrow = n, ncol = p)
  X = normal_sample%*%U
  return(X)
}
sample_variational_CS_multiple_m_multi_d_random_design = function(gamma_, m_s, x0_, n, cor, alpha_, cn_ = 1, kern = mat_kern_1_2){
  d = length(x0_)
  if(is.na(cor)){
    print('Uniform Design')
    xns = matrix(runif(n*d, -0.5, 0.5) , nrow = n, ncol = d)
  }else{
    print('Gaussian Design')
    xns = sample_correlated_inputs(n, p = d, rho = cor)  
  }
  
  f0 = f0_xs_multi_d(Xs = xns, alpha = alpha_)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  K_nn = apply(xns, 1, function(x) apply(xns, 1, function(y) kern(x, y, cn = cn)))
  svd_ = svd(K_nn)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  cat('Estimated sigma: ', sigma_est, '\n')
  
  vp_fit = vp_mean_variance_multiple_m_multi_d(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                               y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  vp_CS = cbind(vp_mean - sqrt(vp_var)*z_gamma, vp_mean+sqrt(vp_var)*z_gamma)
  return(vp_CS)
}

compute_NLPD_multi_d_random_design = function(gamma_, m_s, x0_, n, cor, alpha_, cn_ = 1, kern = mat_kern_1_2){
  d = length(x0_)
  if(is.na(cor)){
    print('Uniform Design')
    xns = matrix(runif(n*d, -0.5, 0.5) , nrow = n, ncol = d)
  }else{
    print('Gaussian Design')
    xns = sample_correlated_inputs(n, p = d, rho = cor)  
  }
  
  f0 = f0_xs_multi_d(Xs = xns, alpha = alpha_)
  y = f0 + rnorm(n, mean = 0, sd = 1)
  K_nn = apply(xns, 1, function(x) apply(xns, 1, function(y) kern(x, y, cn = cn)))
  svd_ = svd(K_nn)
  sigma_est = maximise_sigma_known_svd(y, svd_)
  cat('Estimated sigma: ', sigma_est, '\n')
  
  vp_fit = vp_mean_variance_multiple_m_multi_d(x_=x0_, m_s=m_s, svd_=svd_, k_=kern, xns_=xns,
                                               y_=y, cn_ = cn_, sigma = sigma_est)
  vp_mean = vp_fit[,1]
  vp_var = vp_fit[,2]
  
  z_gamma = qnorm((1+gamma_)/2)
  NLPD = -log(dnorm(0, vp_mean, sd = sqrt(vp_var)))
  return(NLPD)
}



estimate_cov_len_bias = function(gamma_, m_, x0_, svd_, alpha_, num_replicates = 1000, cn_ = 1, rand_design = FALSE){
  cat('m_: ',m_,'\n')
  K = 401
  f0_ks = (1:K)^(-1/2 - alpha_)
  f0_x0 = as.numeric(compute_fn_from_coefs_x(f0_ks, xs = x0))
  test_intervals = t(replicate(num_replicates,
                               sample_variational_CS(gamma_, m_ = m_,
                                                     x0_ = x0_, svd_ = svd_,
                                                     alpha_ = alpha_, cn_ = cn_)))
  
  cov = mean(test_intervals[,1] <= f0_x0 & test_intervals[,2] >= f0_x0)
  len = mean(test_intervals[,2] - test_intervals[,1])
  bias = mean((test_intervals[,1]+test_intervals[,2])/2 - f0_x0)
  
  return(c(cov, len, bias))
}

estimate_cov_len_bias_multiple_m = function(gamma_, m_s, x0_, svd_, alpha_, num_replicates = 1000, cn_ = 1,
                                            kern = bm_kern, xns, noise = 'gaussian'){
  n = dim(svd_$u)[1]
  # K = 401
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0_x0 = as.numeric(compute_fn_from_coefs_x(f0_ks, xs = x0_))
  f0_x0 = f0_xs(xs = x0_, alpha = alpha_)
  xns = c(1:n)/(n+1/2)
  test_intervals = replicate(num_replicates,
                             sample_variational_CS_multiple_m(gamma_, m_s = m_s,
                                                              x0_ = x0_, svd_ = svd_,
                                                              alpha_ = alpha_, cn_ = cn_, kern = kern,
                                                              xns = xns, noise = noise))
  NLPDs = replicate(num_replicates,
                    compute_NLPD(gamma_, m_s = m_s,
                                 x0_ = x0_, svd_ = svd_,
                                 alpha_ = alpha_, cn_ = cn_, kern = kern,
                                 xns = xns, noise = noise))
  
  cat('alpha: ', alpha_, '\n')
  cat('f_0(x_0): ', f0_x0, '\n')
  ret_mat = matrix(nrow = length(m_s), ncol = 5)
  for(ix in 1:length(m_s)){
    m_intervals = t(test_intervals[ix,,])
    ret_mat[ix,1] = mean(m_intervals[,1] <= f0_x0 & m_intervals[,2] >= f0_x0)
    ret_mat[ix,2] = mean(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,3] = sqrt(mean(((m_intervals[,1]+m_intervals[,2])/2 - f0_x0)^2))
    ret_mat[ix, 4] = mean(NLPDs[ix,])
    ret_mat[ix, 5] = sd(NLPDs[ix,])
  }
  return(ret_mat)
}

estimate_cov_len_bias_multiple_m_random_design = function(gamma_, m_s, x0_, n, alpha_, num_replicates = 1000, cn_ = 1,
                                                          kern = bm_kern, noise = 'gaussian'){
  f0_x0 = f0_xs(xs = x0_, alpha = alpha_)
  test_intervals = replicate(num_replicates,
                             sample_variational_CS_multiple_m_random_design(gamma_, m_s = m_s,
                                                                            x0_ = x0_, n = n,
                                                                            alpha_ = alpha_, cn_ = cn_, kern = kern, noise = noise))
  NLPDs = replicate(num_replicates,
                    compute_NLPD_random_design(gamma_, m_s = m_s,
                                               x0_ = x0_, n = n,
                                               alpha_ = alpha_, cn_ = cn_, kern = kern, noise = noise))
  
  cat('alpha: ', alpha_, '\n')
  cat('f_0(x_0): ', f0_x0, '\n')
  ret_mat = matrix(nrow = length(m_s), ncol = 6)
  for(ix in 1:length(m_s)){
    m_intervals = t(test_intervals[ix,,])
    ret_mat[ix,1] = mean(m_intervals[,1] <= f0_x0 & m_intervals[,2] >= f0_x0)
    ret_mat[ix,2] = mean(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,3] = sd(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,4] = sqrt(mean(((m_intervals[,1]+m_intervals[,2])/2 - f0_x0)^2))
    ret_mat[ix,5] = mean(NLPDs[ix,])
    ret_mat[ix,6] = sd(NLPDs[ix,])
  }
  return(ret_mat)
}

estimate_cov_len_bias_multiple_m_multi_d = function(gamma_, m_s, x0_, svd_, alpha_, num_replicates = 1000, cn_ = 1,
                                                    kern = mat_kern_1_2, xns){
  n = dim(svd_$u)[1]
  # K = 401
  # f0_ks = (1:K)^(-1/2 - alpha_)
  # f0_x0 = sum(sapply(x0_, function(x) compute_fn_from_coefs_x(f0_ks, xs = x)))
  f0_x0 = f0_xs_multi_d(Xs = x0_, alpha = alpha_)
  
  # TODO: ESTIMATE CN HERE
  
  test_intervals = replicate(num_replicates,
                             sample_variational_CS_multiple_m_multi_d(gamma_, m_s = m_s,
                                                                      x0_ = x0_, svd_ = svd_,
                                                                      alpha_ = alpha_, cn_ = cn_, kern = kern,
                                                                      xns = xns))
  
  NLPDs = replicate(num_replicates,
                    compute_NLPD_multi_d(gamma_, m_s = m_s,
                                         x0_ = x0_, svd_ = svd_,
                                         alpha_ = alpha_, cn_ = cn_, kern = kern,
                                         xns = xns))
  
  ret_mat = matrix(nrow = length(m_s), ncol = 6)
  for(ix in 1:length(m_s)){
    m_intervals = t(test_intervals[ix,,])
    ret_mat[ix,1] = mean(m_intervals[,1] <= f0_x0 & m_intervals[,2] >= f0_x0)
    ret_mat[ix,2] = mean(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,3] = sd(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,4] = sqrt(mean(((m_intervals[,1]+m_intervals[,2])/2 - f0_x0)^2))
    ret_mat[ix,5] = mean(NLPDs[ix,])
    ret_mat[ix,6] = sd(NLPDs[ix,])
  }
  return(ret_mat)
}
estimate_cov_len_bias_multiple_m_multi_d_emp_bayes = function(gamma_, m_s, x0_, svd_, alpha_, num_replicates = 1000, cn_ = 1,
                                                              kern = mat_kern_1_2, xns){
  n = dim(svd_$u)[1]
  f0_x0 = f0_xs_multi_d(Xs = x0_, alpha = alpha_)
  y = f0_x0 + rnorm(n)
  print('Fitting cn')
  cn = maximise_cn(y, xns, kern)
  cat('maximum cn: ', cn, 'gamma: ', -(1 + log(cn, n))/(2*log(cn, n)))
  gamma_hat = -(1 + log(cn, n))/(2*log(cn, n))
  # cn = n^{gamma_hat - 1}
  
  
  K_nn = apply(xns, 1, function(x) apply(xns, 1, function(y) kern(x, y, cn = cn)))
  svd = svd(K_nn)
  
  test_intervals = replicate(num_replicates,
                             sample_variational_CS_multiple_m_multi_d(gamma_, m_s = m_s,
                                                                      x0_ = x0_, svd_ = svd,
                                                                      alpha_ = alpha_, cn_ = cn, kern = kern,
                                                                      xns = xns))
  NLPDs = replicate(num_replicates,
                    compute_NLPD_multi_d(gamma_, m_s = m_s,
                                         x0_ = x0_, svd_ = svd_,
                                         alpha_ = alpha_, cn_ = cn, kern = kern,
                                         xns = xns))
  
  ret_mat = matrix(nrow = length(m_s), ncol = 5)
  for(ix in 1:length(m_s)){
    m_intervals = t(test_intervals[ix,,])
    ret_mat[ix,1] = mean(m_intervals[,1] <= f0_x0 & m_intervals[,2] >= f0_x0)
    ret_mat[ix,2] = mean(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,3] = sqrt(mean(((m_intervals[,1]+m_intervals[,2])/2 - f0_x0)^2))
    ret_mat[ix,4] = mean(NLPDs[ix,])
    ret_mat[ix,5] = sd(NLPDs[ix,])
  }
  cat('maximum cn: ', cn, 'gamma: ', -(1 + log(cn, n))/(2*log(cn, n)))
  return(ret_mat)
}


estimate_cov_len_bias_multiple_m_multi_d_random_design = function(gamma_, m_s, x0_, n, cor, alpha_, num_replicates = 1000, cn_ = 1,
                                                                  kern = mat_kern_1_2){
  f0_x0 = f0_xs_multi_d(Xs = x0_, alpha = alpha_)
  cat('alpha: ', alpha_, '\n')
  cat('f_0(x_0): ', f0_x0, '\n')
  print(f0_x0)
  test_intervals = replicate(num_replicates,
                             sample_variational_CS_multiple_m_multi_d_random_design(gamma_, m_s = m_s,
                                                                                    x0_ = x0_, n, cor,
                                                                                    alpha_ = alpha_, cn_ = cn_, kern = kern))
  NLPDs = replicate(num_replicates,
                    compute_NLPD_multi_d_random_design(gamma_, m_s = m_s,
                                                       x0_ = x0_, n, cor,
                                                       alpha_ = alpha_, cn_ = cn_, kern = kern))
  
  ret_mat = matrix(nrow = length(m_s), ncol = 6)
  for(ix in 1:length(m_s)){
    m_intervals = t(test_intervals[ix,,])
    ret_mat[ix,1] = mean(m_intervals[,1] <= f0_x0 & m_intervals[,2] >= f0_x0)
    ret_mat[ix,2] = mean(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,3] = sd(m_intervals[,2] - m_intervals[,1])
    ret_mat[ix,4] = sqrt(mean(((m_intervals[,1]+m_intervals[,2])/2 - f0_x0)^2))
    ret_mat[ix,5] = mean(NLPDs[ix,])
    ret_mat[ix,6] = sd(NLPDs[ix,])
  }
  return(ret_mat)
}

r_m = function(x_, xns_, m_, svd_, k_){
  eta_ks = 1/(svd_$d[1:m_] + 1)
  V_m = svd_$u[,1:m_]
  A = V_m %*% diag(eta_ks) %*% t(V_m)
  k_n_x = sapply(xns_, function(y) k_(y, x_))
  return(A%*%k_n_x)
}

#### LIKELIHOOD FUNCTIONS
lmlikelihood_known_svd = function(y, sigma, svd_){
  n = length(y)
  K_n_inv = svd_$v %*% diag(1/(sigma^2 + svd_$d)) %*% t(svd_$v)
  return( -(1/2)*t(y) %*% K_n_inv %*% y -(1/2)*sum(log(svd_$d + sigma^2)))
}
maximise_sigma_known_svd = function(y, svd_){
  sigmas = seq(0.1, 2, by = 0.1)
  results = sapply(sigmas, function(s) lmlikelihood_known_svd(y, s, svd_))
  max_sigma = sigmas[which.max(results)]
  return(max_sigma)
}

lmlikelihood_unknown_svd = function(y, xns, kern, cn, sigma = 1){
  n = length(y)
  K_nn = apply(xns, 1, function(x) apply(xns, 1, function(y) kern(x, y, cn)))
  svd_ = svd(K_nn)
  K_n_inv = svd_$v %*% diag(1/(sigma^2 + svd_$d)) %*% t(svd_$v)
  return( -(1/2)*t(y) %*% K_n_inv %*% y -(1/2)*sum(log(svd_$d + sigma^2)))
}

maximise_cn = function(y, xns, kern = mat_kern_1_2){
  n = length(y)
  gammas = seq(0.5, 2.0, by = 0.1)
  cns = n^{-1/(1+2*gammas)}
  results = sapply(cns, function(cn) lmlikelihood_unknown_svd(y, xns, kern, cn))
  max_cn = cns[which.max(results)]
  return(max_cn)
}


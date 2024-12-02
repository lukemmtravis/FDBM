rm(list = ls())
setwd('/Users/lmt15/Documents/phd/Variational Inference/London/codes')
source('SVGP_fns.R')


#### Posterior Comparison ####
compute_full_results = function(n){
  kern = bm_kern
  x_obs = seq(0, 1, length.out = n)
  f0 = abs(x_obs - 0.5)
  y = f0 + rnorm(length(x_obs), sd = 0.1)
  gamma = 0.5
  cn = (n+1/2)^{(1-2*gamma)/(1+2*gamma)}
  
  xs = seq(0.2, 0.8, by = 0.01)
  f0_xs = abs(xs - 0.5)
  K_nn = sapply(x_obs, function(x) sapply(x_obs, function(y) bm_kern(x, y, cn)))
  svd_decomp = svd(K_nn)
  
  m_bad = n^{2/8}
  m_good = n^{3/4}
  vp_mean_variance(0.1, n, svd_decomp, bm_kern, x_obs, y, cn)
  full_post = t(sapply(xs, function(x) vp_mean_variance(x, n, svd_decomp, bm_kern, x_obs, y, cn)))
  bad_var_post = t(sapply(xs, function(x) vp_mean_variance(x, m_bad, svd_decomp, bm_kern, x_obs, y, cn)))
  good_var_post = t(sapply(xs, function(x) vp_mean_variance(x, m_good, svd_decomp, bm_kern, x_obs, y, cn)))
  
  full_post = data.frame(full_post)
  colnames(full_post) = c('Mean', 'Variance')
  full_post['Posterior'] = 'CFull Posterior'
  full_post['x'] = xs
  
  bad_var_post = data.frame(bad_var_post)
  colnames(bad_var_post) = c('Mean', 'Variance')
  bad_var_post['Posterior'] = 'Am = 5'
  # bad_var_post['Posterior'] = 'ASGPR (m = n^{1/4})'
  bad_var_post['x'] = xs
  
  good_var_post = data.frame(good_var_post)
  colnames(good_var_post) = c('Mean', 'Variance')
  good_var_post['Posterior'] = 'Bm = 106'
  good_var_post['x'] = xs
  
  
  # good_var_post['Mean'] = good_var_post['Mean'] + 0.005
  full_results = rbind(full_post, bad_var_post, good_var_post)
  return(list(full_results = full_results, x_obs = x_obs, f0 = f0, y = y,
              xs = xs, f0_xs = f0_xs))
}
set.seed(1)
n=1000
temp = compute_full_results(n)
full_results = temp$full_results
x_obs = temp$x_obs
y = temp$y
f0 = temp$f0_xs
xs = temp$xs

ggplot(data = full_results) +
  geom_line(data = full_results, aes(x = x, y = Mean, color = Posterior), alpha = 1) +
  geom_ribbon(data= full_results, aes(x = x, ymin = Mean - 1.96*Variance, ymax = Mean + 1.96*Variance, fill = Posterior), color = NA, alpha = 0.1) +
  geom_line(data = data.frame(xs, f0), aes(x = xs, y = f0), color = 'black', alpha = 0.75) +
  xlim(0.2, 0.8) +
  # facet_wrap(~Posterior)
  # ggtitle('Comparison of the Posteriors') +
  ylim(-0.1, 0.45) +
  scale_color_discrete(labels=c('SGPR (m = 5)', 'SGPR (m = 106)', 'GP')) +
  scale_fill_discrete(labels=c('SGPR (m = 5)', 'SGPR (m = 106)', 'GP')) +
  theme(legend.position="bottom", strip.background = element_blank(), strip.text.x = element_blank())

#### Investigate kernels from eigenvectors ####
n = 100
# xns = 5*((c(1:n)  - n/2)/(n + 0.5))
xns = c(1:n)/(n+0.5)
# xns = sort(runif(n))
# xns = c(1:n)/(n+0.5)
start_ix = 1
end_ix = 6

test_kern = bm_kern
K = sapply(xns, function(x) sapply(xns, function(y) test_kern(x, y)) )
svd_decomp = svd(K)
vecs = data.frame(coord = c(1:n), kern = 'bm_dat', svd_decomp$u[,(start_ix:end_ix)])
colnames(vecs) = c(c('coord', 'Kernel'), c(start_ix:end_ix))
bm_dat = gather(vecs, key = 'vector', value = 'v_ij', -c(coord, Kernel))

test_kern = mat_kern_3_2
K = sapply(xns, function(x) sapply(xns, function(y) test_kern(x, y)) )
svd_decomp = svd(K)
vecs = data.frame(coord = c(1:n), kern = 'mat_3_2', svd_decomp$u[,(start_ix:end_ix)])
colnames(vecs) = c(c('coord', 'Kernel'), c(start_ix:end_ix))
mat_3_2_dat = gather(vecs, key = 'vector', value = 'v_ij', -c(coord, Kernel))

test_kern = mat_kern_5_2
K = sapply(xns, function(x) sapply(xns, function(y) test_kern(x, y)) )
svd_decomp = svd(K)
vecs = data.frame(coord = c(1:n), kern = 'mat_5_2', svd_decomp$u[,(start_ix:end_ix)])
colnames(vecs) = c(c('coord', 'Kernel'), c(start_ix:end_ix))
mat_5_2_dat = gather(vecs, key = 'vector', value = 'v_ij', -c(coord, Kernel))

test_kern = se_kern
K = sapply(xns, function(x) sapply(xns, function(y) test_kern(x, y)) )
svd_decomp = svd(K)
vecs = data.frame(coord = c(1:n), kern = 'se', svd_decomp$u[,(start_ix:end_ix)])
colnames(vecs) = c(c('coord', 'Kernel'), c(start_ix:end_ix))
se_dat = gather(vecs, key = 'vector', value = 'v_ij', -c(coord, Kernel))

test_kern = mat_kern_1_2
K = sapply(xns, function(x) sapply(xns, function(y) test_kern(x, y)) )
svd_decomp = svd(K)
vecs = data.frame(coord = c(1:n), kern = 'lin', svd_decomp$u[,(start_ix:end_ix)])
colnames(vecs) = c(c('coord', 'Kernel'), c(start_ix:end_ix))
mat_1_2_dat = gather(vecs, key = 'vector', value = 'v_ij', -c(coord, Kernel))

plot_dat = rbind(bm_dat, mat_1_2_dat, mat_3_2_dat, mat_5_2_dat)


ggplot(plot_dat, aes(x = coord, y = -v_ij, col = Kernel)) + geom_point(alpha=0.5) + facet_wrap(~vector, nrow = 2, ncol = 3) +
  geom_hline(yintercept=2/sqrt(2*n + 1), col='red', alpha=0.5) +
  geom_hline(yintercept=-2/sqrt(2*n + 1), col='red', alpha=0.5) +
  scale_color_discrete(name='Kernel', labels = c('rBM', 'Mat-1/2', 'Mat-3/2', 'Mat-5/2')) +
  ylab(TeX(r'($v_{k}^l$)')) +
  xlab(TeX(r'($l$)'))

#### Fast matrix inversion routine ####
evec_from_psi = function(psi, n){
  vec = (2/sqrt(2*n+1))*sin(c(1:n)*psi)
}

manual_svd = function(n, cn){
  N = (n+1/2)/cn
  psis = (c(1:n) - 1/2)/(n+1/2) * pi
  mus = 1/(2*N*(1-cos(psis)))
  vecs = sapply(psis, function(p) evec_from_psi(p, n) )
  # return(list(d = mus, u = vecs))
  return(1)
}

time_est = function(n, num_replicates){
  xns = c(1:n)/(n+0.5)
  kern = bm_kern
  K_nn = sapply(xns, function(x) sapply(xns, function(y) kern(x,y)))
  
  t1_t = Sys.time()
  manual_svd(n, 1) 
  svd(K_nn) 
  t2_t = Sys.time()
  cat('Estimated runtime: ', difftime(t2_t, t1_t, units = 'secs')*num_replicates, '\n')
  
  man_t1 =Sys.time()
  for(i in c(1:num_replicates)){
    manual_svd(n, 1)  
  }
  man_t2 = Sys.time()
  
  t1 =Sys.time()
  for(i in c(1:num_replicates)){
    svd(K_nn)  
  }
  t2 = Sys.time()
  return(c(difftime(man_t2, man_t1, units = 'secs')/num_replicates, 
           difftime(t2, t1, units = 'secs')/num_replicates))
}

num_replicates = 10
ns = c(10,50,100,200,400,800,1600)
results = t(sapply(ns, function(n) time_est(n, num_replicates)))
# results[,1] = results[,1]/results[1,1]
# results[,2] = results[,2]/results[1,2]

plot_dat = data.frame(results, n = ns)
colnames(plot_dat) = c('Direct', 'Traditional', 'n')
plot_dat = plot_dat %>% gather(key = 'method', value = 'time', -n)
ggplot(plot_dat, aes(x = n, y = time, color = method)) +
  geom_point() + 
  geom_line() + 
  scale_color_discrete(name='Method', labels = c('Analytic Expressions', 'Matrix Inversion Routine')) +
  xlab(TeX(r'($ n$)')) +
  ylab('Time (seconds)')
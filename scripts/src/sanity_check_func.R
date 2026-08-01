library(tidyverse)
library(mvtnorm)

calc_loglik_er = function(sim_stats, target_mu, target_sigma) {
  log_lik = dnorm(x = sim_stats["Mean Degree"], mean = target_mu, sd = sqrt(target_sigma), log = TRUE)
  return(log_lik)
}

# =========================================================
# Approach 1 (BM style)
# =========================================================
generate_network_bm_er = function(params, static_data) {
  alpha = unname(params['alpha'])
  
  p = plogis(alpha)
  n_pairs = length(static_data$rows)
  
  edge_vec = rbinom(n_pairs, 1, p)
  connected_idx = which(edge_vec == 1)
  
  if (length(connected_idx) > 0) {
    return(list(
      rows = static_data$rows[connected_idx],
      cols = static_data$cols[connected_idx],
      has_edge = TRUE
    ))
  } else {
    return(list(has_edge = FALSE))
  }
}

calc_sim_stats_bm_er = function(network, static_data) {
  if (!network$has_edge) return(c("Mean Degree" = 0))
  
  all_indices = c(network$rows, network$cols)
  degree_counts = tabulate(all_indices, nbins = static_data$N)
  
  return(c("Mean Degree" = mean(degree_counts)))
}

# =========================================================
# Approach 2 (M2 style)
# =========================================================
calc_prob_m2_er = function(params, static_data) {
  alpha = unname(params['alpha'])
  p = plogis(alpha)
  n_pairs = length(static_data$rows)
  
  prob_vec = rep(p, n_pairs)
  return(prob_vec)
}

calc_expected_stats_m2_er = function(prob_vec, static_data) {
  N = static_data$N
  p = prob_vec[1] 
  
  expected_mean_degree = (N - 1) * p
  
  return(c("Mean Degree" = expected_mean_degree))
}
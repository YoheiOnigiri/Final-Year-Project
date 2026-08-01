library(tidyverse)
library(mvtnorm)

run_mcmc = function(n_iter, current_params, proposal_cov,
                    static_data, target_mu, target_sigma,
                    generate_func, calc_stats_func, loglik_func) {
  
  # setup
  n_params = length(current_params)
  chain = matrix(NA, nrow = n_iter, ncol = n_params)
  colnames(chain) = names(current_params)
  
  chain[1, ] = current_params
  accept_count = 0
  
  init_params = current_params
  
  # the first likelihood
  current_net = generate_func(current_params, static_data)
  current_stats = calc_stats_func(current_net, static_data)
  current_loglik = loglik_func(current_stats, target_mu, target_sigma)
  
  for (i in 2 : n_iter) {
    
    # A. Propose
    # Metropolis Hastings
    # Sample the step from multivariate normal distribution (dimension = 5)
    proposed_step = as.vector(rmvnorm(n = 1, mean = rep(0, n_params), sigma = proposal_cov))
    proposed_params = current_params + proposed_step
    
    # handle edge cases
    # Sigma cannot be negative
    
    # if (proposed_params['sigma'] <= 0) {
    if (proposed_params['sigma_re'] <= 0) {
      proposed_loglik = -Inf
    }
    else {
      
      # B. Generate the Network
      sim_net_new = generate_func(proposed_params, static_data)
      
      # C. Evaluate the Simulated Network
      sim_stats_new = calc_stats_func(sim_net_new, static_data)
      proposed_loglik = loglik_func(sim_stats_new, target_mu, target_sigma)
      
    }
    
    # D. Accept / Reject
    # Compute the ratio of two log likelihood -> log(L_new / L_old)
    log_ratio = proposed_loglik - current_loglik
    
    # if L_new > L_old -> always accept
    # if L_new < L_old -> accpet stochastically according to the ratio between the two
    
    # Compute the acceptance rate
    if (is.finite(proposed_loglik) && log(runif(1)) < log_ratio) {
      # Accept!
      current_params = proposed_params
      current_loglik = proposed_loglik
      accept_count = accept_count + 1
    } else {
      # Reject
      # No change
    }
    
    # E. Save the result
    chain[i, ] = current_params
    
    # # Show progress
    # if (i %% 1000 == 0) {
    #   cat('Iteration:', i,'\n', ' Current p:', round(current_params, 2), '\n', 'Acceptance:', accept_count, '\n')
    # }
    
    
  }
  
  return(list(
    # init_params = init_params,
    final_params = current_params,
    # proposal_cov = proposal_cov,
    chain = chain,
    accept_count = accept_count,
    acceptance_rate = accept_count / n_iter,
    n_iter = n_iter
    # target_sigma = target_sigma
  ))
}






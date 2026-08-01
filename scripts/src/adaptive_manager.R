
run_adaptive_mcmc = function(config, model_funcs, static_data, output_dir = NULL) {
  
  # 1. Setup
  current_params = config$init_params
  current_cov = config$init_cov
  target_mu = config$target_mu
  target_sigma = config$target_sigma
  output_dir = config$output_dir
  
  results_history = list()
  
  # 2. Adaptive Loop
  for (i in seq_along(config$inflation_schedule)) {
    
    inflation_scale = config$inflation_schedule[i]
    current_scaling_factor = config$scaling_factors[i]
    
    tic()
    # A. MCMC
    mcmc_result = run_mcmc(
      n_iter = config$n_iter,
      current_params = current_params,
      proposal_cov = current_cov,
      
      generate_func = model_funcs$generate,
      calc_stats_func = model_funcs$calc_stats,
      loglik_func = model_funcs$loglik,
      
      static_data = static_data,
      target_mu = target_mu,
      target_sigma = target_sigma * inflation_scale
    )
    stage_time = toc()
    
    # B. Save Log
    stage_config = config
    stage_config$stage = i
    stage_config$inflation_scale = inflation_scale
    
    is_last_stage = (i == length(config$inflation_schedule))
    
    save_mcmc_result(
      result = mcmc_result,
      config = stage_config,
      time_log = stage_time,
      output_dir = output_dir,
      write_md = is_last_stage
    )
    
    # C. Store Result
    results_history[[i]] = list(
      stage = i,
      inflation = inflation_scale,
      chain = mcmc_result$chain,
      acceptance_rate = mcmc_result$acceptance_rate,
      final_params = mcmc_result$final_params
    )
    
    # D. Update Covariance & Params for next iteration
    chain_df = as.data.frame(mcmc_result$chain)
    burn_in_ind   = floor(config$n_iter * config$burn_in)
    valid_chain   = chain_df[-(1 : burn_in_ind), ]
    posterior_cov = cov(valid_chain)
    posterior_cov = posterior_cov + diag(1e-12, ncol(posterior_cov))
    
    current_cov   = posterior_cov * current_scaling_factor * config$safety_margin
    current_params = mcmc_result$final_params
    
    config$init_params = current_params
    config$init_cov = current_cov
    
    cat(sprintf("Stage %d (Inflation: %.1f) Acceptance Rate: %.2f%%\n", 
                i, inflation_scale, mcmc_result$acceptance_rate * 100))
  }
  
  return(list(
    history = results_history,
    final_params = current_params,
    final_cov = current_cov
  ))
}













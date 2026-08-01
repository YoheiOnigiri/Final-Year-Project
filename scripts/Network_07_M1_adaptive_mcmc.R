# 1. Setup =====================================================================

library(tidyverse)
library(here)
library(tictoc)
library(mvtnorm)

func_path = here('scripts', 'src')
rds_path = here('data', 'processed')

# General functions
source(here(func_path, 'mcmc_core.R'))
source(here(func_path, 'adaptive_manager.R'))
source(here(func_path, 'log_utils.R'))
# For this model
source(here(func_path, 'model_BM.R'))

# Necessary Data
target_stats = readRDS(here(rds_path, 'BM_target_stats.rds'))
nodes = readRDS(here(rds_path, 'simulation_nodes.rds'))

# Static Data for Network Generation
static_data = get_static_data_bm(nodes)


# 2. Hyper Parameters ==========================================================

# --- MCMC Basic Settings ---
n_iter        = 200000
burn_in_ratio = 0.2
safety_margin = 0.75
init_const    = 0.01 ^ 2

# --- Inflation & Scaling Schedule ---
optimal_scale      = (2.38 ^ 2) / 5
inflation_schedule = c(50, 20, 10, 5, 2, 1)
next_schedule      = c(inflation_schedule[-1], NA)
shrinkage_ratios   = next_schedule / inflation_schedule

scaling_factors    = optimal_scale * shrinkage_ratios
scaling_factors[is.na(scaling_factors)] = 1.0

# --- Initial Covariance Matrix ---
init_proposal_sd = c(
  alpha  = 0.07,     
  theta  = 0.002,     
  phi_mm = 0.08,    
  phi_ff = 0.08,    
  sigma  = 0.03
)

init_proposal_cov = diag(init_proposal_sd ^ 2) * init_const 
rownames(init_proposal_cov) = names(init_proposal_sd)
colnames(init_proposal_cov) = names(init_proposal_sd)

# --- Initial Params ---
init_params = c(
  alpha  = -4.1,  
  theta  = -0.040,
  phi_mm = 0.25,
  phi_ff = 0.29,
  sigma  = 0.68
)

# 3. Config Setup ==============================================================

target_mu = target_stats$mu
target_sigma = target_stats$sigma

config = list(
  theme = 'BM',
  name = 'test3',
  note = 'Test run for new bootstrap',
  
  n_iter = n_iter,
  burn_in = burn_in_ratio,
  
  inflation_schedule = inflation_schedule,
  scaling_factors = scaling_factors,
  
  safety_margin = safety_margin,
  
  init_params = init_params,
  init_cov = init_proposal_cov,
  
  target_mu = target_mu,
  target_sigma = target_sigma
)

# 4. Model Functions Packing ===================================================

model_funcs = list(
  generate = generate_network_bm,
  calc_stats = calc_sim_stats_bm,
  loglik = calc_loglik_bm
)

# 5. Run Adaptive MCMC =========================================================

results = run_adaptive_mcmc(
  config = config,
  model_funcs = model_funcs,
  static_data = static_data,
  output_dir = rds_path
)

final_chain = results$history[[length(results$history)]]$chain
par(mfrow=c(2,3))
for(p in colnames(final_chain)){
  plot(final_chain[,p], type='l', main=p, ylab="Value", xlab="Iter")
}

# 6. Check =====================================================================
set.seed(100)
sim_net_test = model_funcs$generate(config$init_params, static_data)
sim_stats_test = model_funcs$calc_stats(sim_net_test, static_data)

stats_names = c("Mean Degree", "Log Variance", "Mean Age Diff", "Prop MM", "Prop FF")

comparison = data.frame(
  Stats_Name = stats_names,
  Generated  = as.numeric(sim_stats_test),
  Target_Mu  = as.numeric(config$target_mu),
  Target_SD  = sqrt(diag(config$target_sigma)),
  Diff_Sigma = abs(as.numeric(sim_stats_test) - as.numeric(config$target_mu)) / sqrt(diag(config$target_sigma))
)

print(comparison)

ll = model_funcs$loglik(sim_stats_test, config$target_mu, config$target_sigma)
cat("\nLog Likelihood (Inflation 1.0):", ll, "\n")
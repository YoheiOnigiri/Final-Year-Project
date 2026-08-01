library(tidyverse)
library(mvtnorm)

# 1. Static Data ---------------------------------------------------------------
get_static_data_bm = function(nodes) {
  N = nrow(nodes)
  
  # Get the indices of upper triangle matrix
  upper_tri_ind = which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
  
  # vectorize the necessary information
  ages = nodes$part_age
  genders = nodes$part_gender
  is_male = (genders == 'M')
  is_female = (genders == 'F')
  
  rows = upper_tri_ind[, 1]
  cols = upper_tri_ind[, 2]
  
  return(list(
    N = N,
    rows = rows,
    cols = cols,
    # Age Difference
    age_diff = abs(ages[rows] - ages[cols]),
    # Dummy MM
    male_male = (is_male[rows] & is_male[cols]) * 1,
    # Dummy FF
    female_female = (is_female[rows] & is_female[cols]) * 1,
    # For summary statistics
    part_ages = ages,
    part_genders = genders
  ))
}

# 2. Network Generation --------------------------------------------------------

generate_network_bm = function(params, static_data) {
  
  # ----- parameters ---------------------
  
  alpha = params['alpha']
  theta = params['theta']
  phi_mm = params['phi_mm']
  phi_ff = params['phi_ff']
  
  # sigma = params['sigma']
  sigma_re = params['sigma_re']
  
  N = static_data$N
  
  # ----- Random Effect -------------------
  
  # epsilon = rnorm(N, mean = 0, sd = sigma)
  epsilon = rnorm(N, mean = 0, sd = sigma_re)
  
  # calculate the necessary pairs instead of creating an entire matrix
  re_sum_vec = epsilon[static_data$rows] + epsilon[static_data$cols]
  
  # ----- Logit ----------------------------
  
  logit_P = alpha +
    theta * static_data$age_diff +
    phi_mm * static_data$male_male +
    phi_ff * static_data$female_female +
    re_sum_vec
  
  # ----- Logit Space Correction -----------
  logit_P_corrected = logit_P - log(N / 1000)
  
  # ----- Convert Logit into Probability ---
  
  prob_vec = plogis(logit_P_corrected)
  
  # ----- Create Edges ---------------------
  
  edge_vec = rbinom(n = length(prob_vec), size = 1, prob = prob_vec)
  
  # ----- Return indices where edge exists -
  
  # Indices of connected pairs
  connected_indices = which(edge_vec == 1)
  
  if (length(connected_indices) > 0) {
    # only return the indices of nodes that are connected
    return(list(
      rows = static_data$rows[connected_indices],
      cols = static_data$cols[connected_indices],
      has_edge = TRUE
    ))
  } else {
    return(list(has_edge = FALSE))
  }
}

# 3. Calculate Summary Statistics ----------------------------------------------

calc_sim_stats_bm = function(network, static_data) {
  # network -> edge_list
  
  # Handle edge case
  if (!network$has_edge) {
    return(c(
      "Mean Degree" = 0,
      "Log Variance" = -Inf, # log(-Inf) = 0
      "Mean Age Diff" = 100, # give large number as a penalty
      "Prop MM" = 0,
      "Prop FF" = 0
    ))
  }
  
  # 1. Degree and Variance
  # Calculate the degree of all nodes
  all_indices = c(network$rows, network$cols)
  degree_counts = tabulate(all_indices, nbins = static_data$N)
  
  mean_degree = mean(degree_counts)
  var_degree = var(degree_counts)
  if (var_degree == 0) {
    var_degree = 1e-6
  }
  
  # 2. Joining the information of the nodes & Calculate summary statistics
  # Use index reference instead of dplyr to make it faster
  
  # The indices of connected pairs
  p_rows = network$rows
  p_cols = network$cols
  
  # Age Difference
  ages = static_data$part_ages
  age_diffs = abs(ages[p_rows] - ages[p_cols])
  mean_age_diff = mean(age_diffs)
  
  # Gender Homophily
  genders = static_data$part_genders
  gens_row = genders[p_rows]
  gens_col = genders[p_cols]
  
  # Prop MM
  prop_mm = mean(gens_row == 'M' & gens_col == 'M')
  # Prop FF
  prop_ff = mean(gens_row == 'F' & gens_col == 'F')
  
  return(c(
    mean_degree,
    log(var_degree),
    mean_age_diff,
    prop_mm,
    prop_ff
  ))
}

# 4. Calculate Log Likelihood --------------------------------------------------

calc_loglik_bm = function(sim_stats, target_mu, target_sigma) {
  log_lik = dmvnorm(x = sim_stats, mean = target_mu, sigma = target_sigma, log = TRUE)
  return(log_lik)
}
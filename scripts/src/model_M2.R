library(tidyverse)
library(mvtnorm)

# 1. Static Data ---------------------------------------------------------------
get_static_data_m2 = function(nodes) {
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
  
  is_school = (ages <= 20 & ages >= 5)
  school_sum = (is_school[rows] & is_school[cols]) * 1
  
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
    # School or not
    school_sum = school_sum,
    # For summary statistics
    part_ages = ages,
    part_genders = genders
  ))
}

# 2. Network Generation --------------------------------------------------------

generate_network_m2 = function(params, static_data) {
  
  # ----- parameters ---------------------
  
  alpha = params['alpha']
  gamma_school = params['gamma_school']
  w = params['w']
  lambda = params['lambda']
  sigma_hhd = params['sigma_hhd']
  mu = params['mu']
  
  phi_mm = params['phi_mm']
  phi_ff = params['phi_ff']
  sigma_re = params['sigma_re']
  
  N = static_data$N
  
  # ----- Random Effect -------------------
  
  epsilon = rnorm(N, mean = 0, sd = sigma_re)
  
  # calculate the necessary pairs instead of creating an entire matrix
  re_sum_vec = epsilon[static_data$rows] + epsilon[static_data$cols]
  
  # ----- Age Difference Distribution -----
  age_diff = static_data$age_diff
  
  # same age (exponential distribution)
  dens_exp = lambda * exp(-lambda * age_diff)
  
  # household (normal distribution)
  dens_norm = (1 / (sqrt(2 * pi) * sigma_hhd)) * exp(-((age_diff - mu)^2) / (2 * sigma_hhd^2))
  
  # adding the two distributions up
  mix_density = w * dens_exp + (1 - w) * dens_norm
  
  # prevent log(0)
  mix_density = pmax(mix_density, 1e-10)
  
  # ----- Logit ----------------------------
  
  logit_P = alpha +
    gamma_school * static_data$school_sum +
    log(mix_density) +
    phi_mm * static_data$male_male +
    phi_ff * static_data$female_female +
    re_sum_vec
  
  # ----- Convert Logit into Probability ---
  
  prob_vec = plogis(logit_P)
  
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

calc_sim_stats_m2 = function(network, static_data) {
  # network -> edge_list
  
  # Handle edge case
  if (!network$has_edge) {
    return(c(
      "Mean Degree" = 0,
      "Mean Degree School" = 0,
      "Log Variance" = -Inf, # log(-Inf) = 0
      "Prop MM" = 0,
      "Prop FF" = 0,
      "Prop_age_0_5" = 0,
      "Prop_age_20_35" = 0
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
  
  # School check
  school_indices = which(static_data$part_ages <= 20 & static_data$part_ages >= 5)
  mean_degree_school = mean(degree_counts[school_indices])
  if (is.na(mean_degree_school)) mean_degree_school = 0

    
  # The indices of connected pairs
  p_rows = network$rows
  p_cols = network$cols
  
  # Age Difference
  ages = static_data$part_ages
  age_diffs = abs(ages[p_rows] - ages[p_cols])
  prop_age_0_5 = mean(age_diffs >= 0 & age_diffs <= 5)
  prop_age_20_35 = mean(age_diffs >= 20 & age_diffs <= 35)
  
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
    mean_degree_school,
    log(var_degree),
    prop_mm,
    prop_ff,
    prop_age_0_5,
    prop_age_20_35
  ))
}

calc_sim_stats_m2_2 = function(network, static_data) {
  # network -> edge_list
  
  # Handle edge case
  if (!network$has_edge) {
    return(rep(0, 6))
  }
  
  # 1. Degree and Variance
  # Calculate the degree of all nodes
  all_indices = c(network$rows, network$cols)
  degree_counts = tabulate(all_indices, nbins = static_data$N)
  
  mean_degree = mean(degree_counts)
  
  # 2. Joining the information of the nodes & Calculate summary statistics
  # Use index reference instead of dplyr to make it faster
  
  # The indices of connected pairs
  p_rows = network$rows
  p_cols = network$cols
  
  # Age Difference
  ages = static_data$part_ages
  age_diffs = abs(ages[p_rows] - ages[p_cols])
  prop_age_0_5 = mean(age_diffs >= 0 & age_diffs <= 5)
  prop_age_20_35 = mean(age_diffs >= 20 & age_diffs <= 35)
  
  # Gender Homophily
  genders = static_data$part_genders
  gens_row = genders[p_rows]
  gens_col = genders[p_cols]
  
  # Prop MM
  prop_mm = mean(gens_row == 'M' & gens_col == 'M')
  # Prop FF
  prop_ff = mean(gens_row == 'F' & gens_col == 'F')
  
  # Prop School
  is_school_edge = (ages[p_rows] >= 5 & ages[p_rows] <= 20) & (ages[p_cols] >= 5 & ages[p_cols] <= 20)
  prop_school = mean(is_school_edge)
  
  return(c(
    log(mean_degree),
    prop_mm,
    prop_ff,
    prop_age_0_5,
    prop_age_20_35,
    prop_school
  ))
}

calc_sim_stats_m2_3 = function(network, static_data) {
  # network -> edge_list
  
  # Handle edge case
  if (!network$has_edge) {
    return(rep(0, 9))
  }
  
  # 1. Degree and Variance
  # Calculate the degree of all nodes
  all_indices = c(network$rows, network$cols)
  degree_counts = tabulate(all_indices, nbins = static_data$N)
  
  mean_degree = mean(degree_counts)
  
  prop_degree_le_2 = mean(degree_counts <= 2)
  prop_degree_ge_25 = mean(degree_counts >= 25)
  
  # 2. Joining the information of the nodes & Calculate summary statistics
  # Use index reference instead of dplyr to make it faster
  
  # The indices of connected pairs
  p_rows = network$rows
  p_cols = network$cols
  
  # Age Difference
  ages = static_data$part_ages
  age_diffs = abs(ages[p_rows] - ages[p_cols])
  prop_age_0_5 = mean(age_diffs >= 0 & age_diffs <= 5)
  prop_age_27_32 = mean(age_diffs >= 27 & age_diffs <= 32)
  prop_age_over_40 = mean(age_diffs > 40)
  
  # Gender Homophily
  genders = static_data$part_genders
  gens_row = genders[p_rows]
  gens_col = genders[p_cols]
  
  # Prop MM
  prop_mm = mean(gens_row == 'M' & gens_col == 'M')
  # Prop FF
  prop_ff = mean(gens_row == 'F' & gens_col == 'F')
  
  # Prop School
  is_school_edge = (ages[p_rows] >= 5 & ages[p_rows] <= 20) & (ages[p_cols] >= 5 & ages[p_cols] <= 20)
  prop_school = mean(is_school_edge)
  
  return(c(
    log(mean_degree),
    prop_degree_le_2,
    prop_degree_ge_25,
    prop_mm,
    prop_ff,
    prop_age_0_5,
    prop_age_28_32,
    prop_age_over_40,
    prop_school
  ))
}

# 4. Calculate Log Likelihood --------------------------------------------------

calc_loglik_m2 = function(sim_stats, target_mu, target_sigma) {
  log_lik = dmvnorm(x = sim_stats, mean = target_mu, sigma = target_sigma, log = TRUE)
  return(log_lik)
}
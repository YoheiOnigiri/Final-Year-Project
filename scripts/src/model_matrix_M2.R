library(tidyverse)
library(mvtnorm)

# 1. Static Data ---------------------------------------------------------------
get_static_data_m2_5 = function(nodes, n_school = c(1, 1, 3, 4)) {
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
  
  age_diff = abs(ages[rows] - ages[cols])
  
  is_school = (ages >= 5 & ages <= 22)
  is_school_edge = is_school[rows] & is_school[cols]
  max_age = pmax(ages[rows], ages[cols])
  
  # 4 stages of school
  school_pri_low  = (is_school_edge & max_age >= 5 & max_age <= 8 & age_diff <= n_school[1]) * 1
  school_pri_high = (is_school_edge & max_age >= 9 & max_age <= 12 & age_diff <= n_school[2]) * 1
  school_sec      = (is_school_edge & max_age >= 13 & max_age <= 18 & age_diff <= n_school[3]) * 1
  school_ter      = (is_school_edge & max_age >= 19 & max_age <= 22 & age_diff <= n_school[4]) * 1
  
  set.seed(100)
  Z_re = as.numeric(scale(rnorm(N)))
  
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
    school_pri_low = school_pri_low,
    school_pri_high = school_pri_high,
    school_sec = school_sec,
    school_ter = school_ter,
    # For summary statistics
    part_ages = ages,
    part_genders = genders,
    # For Random Effect
    Z_re = Z_re
  ))
}

# 2. Calculate Probability (Network Generation)--------------------------------------------------------

calc_prob_m2_3 = function(params, static_data) {
  
  # ----- parameters ---------------------
  
  alpha = params['alpha']
  gamma_pri_low = params['gamma_pri_low']
  gamma_pri_high = params['gamma_pri_high']
  gamma_sec = params['gamma_sec']
  gamma_ter = params['gamma_ter']
  w = params['w']
  lambda = params['lambda']
  sigma_hhd = params['sigma_hhd']
  mu = params['mu']
  
  phi_mm = params['phi_mm']
  phi_ff = params['phi_ff']
  sigma_re = params['sigma_re']
  
  N = static_data$N
  
  # ----- Random Effect -------------------
  
  # epsilon = rnorm(N, mean = 0, sd = sigma_re)
  # 
  # # calculate the necessary pairs instead of creating an entire matrix
  # re_sum_vec = epsilon[static_data$rows] + epsilon[static_data$cols]
  
  # ----- New Random Effect ---------------
  epsilon = static_data$Z_re * sigma_re
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
    gamma_pri_low * static_data$school_pri_low +
    gamma_pri_high * static_data$school_pri_high +
    gamma_sec * static_data$school_sec +
    gamma_ter * static_data$school_ter +
    log(mix_density) +
    phi_mm * static_data$male_male +
    phi_ff * static_data$female_female +
    re_sum_vec
  
  # ----- Convert Logit into Probability ---
  
  prob_vec = plogis(logit_P)
  
  # ----- Return Probability ---------------------
  return(prob_vec)
}

# 3. Calculate Expected Summary Statistics -----------------------------
calc_expected_stats_matrix_m2_6 = function(prob_vec, static_data) {
  
  N = static_data$N
  
  # --- A symmetric NxN Probability Matrix P_mat ---
  P_mat = matrix(0, nrow = N, ncol = N)
  # Fill in both upper and lower triangle
  P_mat[cbind(static_data$rows, static_data$cols)] = prob_vec
  P_mat[cbind(static_data$cols, static_data$rows)] = prob_vec
  
  # --- 1. Degree Stats ---
  
  # Expected Degree of Individual
  d_ind = rowSums(P_mat)
  mean_degree = mean(d_ind)
  
  # Tails (Use Poisson Estimation)
  prop_degree_le_2 = mean(ppois(2, lambda = d_ind))
  prop_degree_ge_25 = mean(1 - ppois(24, lambda = d_ind))
  
  # --- 2. E(Proportions) ---
  
  # Expectation of total contacts (edges) of the network
  total_edges_expected = sum(prob_vec)
  
  # Add up the probabilities that fulfills the condition
  # Prop MM
  prop_mm = sum(prob_vec[static_data$male_male == 1]) / total_edges_expected
  # Prop FF
  prop_ff = sum(prob_vec[static_data$female_female == 1]) / total_edges_expected
  
  # Prop Age 0-5
  age_diffs = static_data$age_diff
  prop_age_0_5 = sum(prob_vec[age_diffs >= 0 & age_diffs <= 5]) / total_edges_expected
  # Prop Age 6-15
  prop_age_6_15 = sum(prob_vec[age_diffs >= 6 & age_diffs <= 15]) / total_edges_expected
  # Prop Age 20-35
  prop_age_20_35 = sum(prob_vec[age_diffs >= 20 & age_diffs <= 35]) / total_edges_expected
  # Prop Age 40-
  prop_age_over_40 = sum(prob_vec[age_diffs > 40]) / total_edges_expected
  
  # Prop School
  prop_school_pri_low = sum(prob_vec[static_data$school_pri_low == 1]) / total_edges_expected
  prop_school_pri_high = sum(prob_vec[static_data$school_pri_high == 1]) / total_edges_expected
  prop_school_sec = sum(prob_vec[static_data$school_sec == 1]) / total_edges_expected
  prop_school_ter = sum(prob_vec[static_data$school_ter == 1]) / total_edges_expected
  
  return(c(
    log(mean_degree),
    prop_degree_le_2,
    prop_degree_ge_25,
    prop_mm,
    prop_ff,
    prop_age_0_5,
    prop_age_6_15,
    prop_age_20_35,
    prop_age_over_40,
    prop_school_pri_low,
    prop_school_pri_high,
    prop_school_sec,
    prop_school_ter
  ))
}

# 4. Calculate Log Likelihood --------------------------------------------------

calc_loglik_m2 = function(sim_stats, target_mu, target_sigma) {
  log_lik = dmvnorm(x = sim_stats, mean = target_mu, sigma = target_sigma, log = TRUE)
  return(log_lik)
}

# 5. Generate Network for Network Simulation -----------------------------------
generate_network_m2_3 = function(params, static_data) {
  
  # ----- parameters ---------------------
  alpha = params['alpha']
  gamma_pri_low = params['gamma_pri_low']
  gamma_pri_high = params['gamma_pri_high']
  gamma_sec = params['gamma_sec']
  gamma_ter = params['gamma_ter']
  w = params['w']
  lambda = params['lambda']
  sigma_hhd = params['sigma_hhd']
  mu = params['mu']
  
  phi_mm = params['phi_mm']
  phi_ff = params['phi_ff']
  sigma_re = params['sigma_re']
  
  N = static_data$N
  
  # ----- New Random Effect ---------------
  epsilon = static_data$Z_re * sigma_re
  re_sum_vec = epsilon[static_data$rows] + epsilon[static_data$cols]
  
  # ----- Age Difference Distribution -----
  age_diff = static_data$age_diff
  
  dens_exp = lambda * exp(-lambda * age_diff)
  dens_norm = (1 / (sqrt(2 * pi) * sigma_hhd)) * exp(-((age_diff - mu)^2) / (2 * sigma_hhd^2))
  
  mix_density = w * dens_exp + (1 - w) * dens_norm
  mix_density = pmax(mix_density, 1e-10)
  
  # ----- Logit ----------------------------
  logit_P = alpha +
    gamma_pri_low * static_data$school_pri_low +
    gamma_pri_high * static_data$school_pri_high +
    gamma_sec * static_data$school_sec +
    gamma_ter * static_data$school_ter +
    log(mix_density) +
    phi_mm * static_data$male_male +
    phi_ff * static_data$female_female +
    re_sum_vec
  
  # ----- Logit Space Correction -----------
  logit_P_corrected = logit_P - log(N / 1000)
  
  # ----- Convert Logit into Probability ---
  prob_vec = plogis(logit_P_corrected)
  
  # ----- Create Edges ---------------------
  edge_vec = rbinom(n = length(prob_vec), size = 1, prob = prob_vec)
  
  connected_indices = which(edge_vec == 1)
  
  if (length(connected_indices) > 0) {
    return(list(
      rows = static_data$rows[connected_indices],
      cols = static_data$cols[connected_indices],
      has_edge = TRUE
    ))
  } else {
    return(list(has_edge = FALSE))
  }
}







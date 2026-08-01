library(tidyverse)
library(here)
library(grid)
library(Matrix)
library(patchwork)
library(deSolve)
library(igraph)

source(here('scripts', 'src', 'model_matrix_M2.R')) 


# run_network_sir updated ----
run_network_sir = function(adj_matrix, num_init_infected = 5, period = 150, rho = 0.1, d = 5, init_method = 'degree') {
  
  N = nrow(adj_matrix)
  node_degrees = Matrix::rowSums(adj_matrix)
  
  # --- Phase 1: Setup ---
  status = rep('S', N)
  days_infected = rep(0, N)
  
  if (init_method == "degree") {
    initial_infected_nodes = sample(1:N, num_init_infected, replace = FALSE, prob = node_degrees)
  } else {
    initial_infected_nodes = sample(1:N, num_init_infected, replace = FALSE)
  }
  
  status[initial_infected_nodes] = 'I'
  
  init_avg_new_contacts = mean(node_degrees[initial_infected_nodes])
  init_avg_current_contacts = init_avg_new_contacts
  
  I_idx = which(status == 'I')
  valid_I_idx = I_idx[node_degrees[I_idx] > 0]
  is_non_S = as.numeric(status != 'S')
  num_non_S_neighbors = as.numeric(adj_matrix %*% is_non_S)
  init_pct_infected = mean(num_non_S_neighbors[valid_I_idx] / node_degrees[valid_I_idx])
  
  results = data.frame(
    time = 0,
    S = sum(status == 'S'),
    I = sum(status == 'I'),
    R = sum(status == 'R'),
    avg_contacts_per_new_case = init_avg_new_contacts,         
    avg_contacts_per_current_case = init_avg_current_contacts, 
    pct_contacts_infected = init_pct_infected                  
  )
  
  # --- Phase 2: Loop ---
  for (day in 1:period) {
    
    # 1. S -> I
    is_infected = as.numeric(status == 'I')
    num_infected_neighbors = as.numeric(adj_matrix %*% is_infected)
    
    susceptible_idx = which(status == 'S')
    prob_infection = 1 - exp(-num_infected_neighbors[susceptible_idx] * rho)
    
    new_infections = rbinom(n = length(susceptible_idx), size = 1, prob = prob_infection)
    new_I_idx = susceptible_idx[new_infections == 1]
    
    # 2. I -> R
    infected_idx = which(status == 'I')
    gamma = 1 / d 
    recoveries = rbinom(n = length(infected_idx), size = 1, prob = gamma)
    new_R_idx = infected_idx[recoveries == 1]
    
    # 3. Update status
    status[new_I_idx] = 'I'
    status[new_R_idx] = 'R'
    
    # --- 4. Record & Calculate Network Metrics ---
    current_I = sum(status == 'I')
    
    if (length(new_I_idx) > 0) {
      avg_new_contacts = mean(node_degrees[new_I_idx])
    } else {
      avg_new_contacts = NA
    }
    
    if (current_I > 0) {
      current_infected_idx = which(status == 'I')
      valid_I_idx = current_infected_idx[node_degrees[current_infected_idx] > 0]
      
      avg_current_contacts = mean(node_degrees[valid_I_idx])
      
      is_non_S = as.numeric(status != 'S')
      num_non_S_neighbors = as.numeric(adj_matrix %*% is_non_S)
      pct_infected = mean(num_non_S_neighbors[valid_I_idx] / node_degrees[valid_I_idx])
    } else {
      avg_current_contacts = NA
      pct_infected = NA
    }
    
    new_data = data.frame(
      time = day,
      S = sum(status == 'S'),
      I = current_I,
      R = sum(status == 'R'),
      avg_contacts_per_new_case = avg_new_contacts,
      avg_contacts_per_current_case = avg_current_contacts,
      pct_contacts_infected = pct_infected
    )
    results = rbind(results, new_data)
    
    if (current_I == 0) break
  }
  
  return(results)
}

# run the simulation (wrapper)
stochastic_sir_simulations = function(adj_list, R0 = 2.5, d = 7, period = 150, 
                                      num_init_infected = 5, seed = 100, 
                                      init_method = 'degree', fixed_rho = NULL) {
  
  num_sims = length(adj_list)
  all_results = list()
  N = nrow(adj_list[[1]])
  set.seed(seed)
  
  for (i in 1 : num_sims) {
    current_adj = adj_list[[i]]
    
    # --- Intervention or not ---
    if (is.null(fixed_rho)) {
      # Baseline : Compute rho from R0
      k_bar = (sum(current_adj)) / N
      current_rho = -log(1 - R0 / k_bar) / d
    } else {
      # Intervention：Use a specific rho
      current_rho = fixed_rho
    }
    # ----------------------
    
    res = run_network_sir(
      adj_matrix = current_adj,
      num_init_infected = num_init_infected,
      period = period,
      rho = current_rho,
      d = d,
      init_method = init_method
    )
    all_results[[i]] = res
  }
  return(all_results)
}

# Model Network generator
generate_sim_networks = function(clean_chain, static_data, N, num_sims = 100, seed = 100,
                                 network_func = generate_network_m2_3) {
  
  sim_adj_list = list()
  set.seed(seed)
  
  sampled_indices = sample(1:nrow(clean_chain), num_sims)
  
  for (i in 1 : num_sims) {
    current_params = clean_chain[sampled_indices[i], ]
    
    sim_network = network_func(current_params, static_data)
    
    if (sim_network$has_edge) {
      adj_matrix = sparseMatrix(
        i = c(sim_network$rows, sim_network$cols),
        j = c(sim_network$cols, sim_network$rows),
        x = 1, dims = c(N, N)
      )
    } else {
      adj_matrix = Matrix(0, nrow = N, ncol = N, sparse = TRUE)
    }
    
    sim_adj_list[[i]] = adj_matrix
  }
  
  cat('Finished!\n')
  return(sim_adj_list)
}

# Random Network (Erdos-Renyi) generator
generate_random_networks = function(sim_adj_list, seed = 100) {
  
  num_sims = length(sim_adj_list)
  N = nrow(sim_adj_list[[1]])
  random_adj_list = list()
  
  upper_tri_indices = which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
  num_possible_edges = nrow(upper_tri_indices)
  
  set.seed(seed)
  
  for (i in 1:num_sims) {

    # 1. i-th model-simulated adjacent matrix
    target_edges = sum(sim_adj_list[[i]]) / 2
    
    # 2. calculate probability p
    p_random = target_edges / num_possible_edges
    
    # 3. roll dice
    er_edges = rbinom(n = num_possible_edges, size = 1, prob = p_random)
    er_connected = upper_tri_indices[er_edges == 1, ]
    
    # 4. create sparse matrix
    adj_random = sparseMatrix(
      i = c(er_connected[, 1], er_connected[, 2]),
      j = c(er_connected[, 2], er_connected[, 1]),
      x = 1, dims = c(N, N)
    )
    
    # contain in list
    random_adj_list[[i]] = adj_random
  }
  
  cat('Finished!')
  return(random_adj_list)
}

# Visualization
plot_single_sir = function(res_list,
                           title = 'Network Stochastic SIR Moddel',
                           subtitle = '100 runs',
                           alpha = 0.1,
                           linewidth = 0.5) {
  
  # 1. Convert the result list to long format
  df_all = bind_rows(res_list, .id = "sim_id") %>%
    pivot_longer(cols = c(S, I, R), names_to = "State", values_to = "Count") %>%
    mutate(State = factor(State, levels = c("S", "I", "R")))
  
  # 2. plot
  p = ggplot(df_all, aes(x = time, y = Count, group = interaction(sim_id, State), color = State)) +
    geom_line(alpha = alpha, linewidth = linewidth) +
    scale_color_manual(values = c("S" = "blue", "I" = "red", "R" = "darkorange")) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time (days)",
      y = "Number of People",
      color = "State"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}

plot_sir_comparison = function(res_model, res_random, 
                               title = "Comparison of Epidemic Dynamics: Model vs Random Network",
                               subtitle = "Both networks have the same average degree") {
  
  # 1. Join the MODEL and RANDOM network simulation results
  df_poly = bind_rows(res_model, .id = "sim_id") %>%
    mutate(Network = "1. MODEL")
  
  df_random = bind_rows(res_random, .id = "sim_id") %>%
    mutate(Network = "2. RANDOM (Erdos-Renyi)")
  
  df_compare = bind_rows(df_poly, df_random) %>%
    pivot_longer(cols = c(S, I, R), names_to = "State", values_to = "Count") %>%
    mutate(State = factor(State, levels = c("S", "I", "R")))
  
  # 2. plot
  p = ggplot(df_compare, aes(x = time, y = Count, group = interaction(sim_id, State), color = State)) +
    geom_line(alpha = 0.1, linewidth = 0.5) +
    scale_color_manual(values = c("S" = "blue", "I" = "red", "R" = "darkorange")) +
    facet_wrap(~ Network) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time (days)",
      y = "Number of People",
      color = "State"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

plot_sir_comparison_vertical = function(res_model, res_random) {
  
  p1 = plot_single_sir(res_model, title = "COVID-19 on MODEL Network", subtitle = "")
  p2 = plot_single_sir(res_random, title = "COVID-19 on RANDOM Network", subtitle = "")

  combined_plot = p1 / p2 + 
    plot_layout(guides = "collect") + 
    plot_annotation(
      title = "Comparison of Epidemic Spread: MODEL vs RANDOM",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  return(combined_plot)
}

plot_sir_overlay = function(res_list, res_det,
                            title = 'Stochastic vs Deterministic SIR',
                            subtitle = 'Overlaid with ODE (thicker lines)',
                            alpha = 0.1,
                            linewidth_stoch = 0.5,
                            linewidth_det = 1.2) {
  
  # 1. Stochastic
  df_stoch = bind_rows(res_list, .id = "sim_id") %>%
    pivot_longer(cols = c(S, I, R), names_to = "State", values_to = "Count") %>%
    mutate(State = factor(State, levels = c("S", "I", "R")))
  
  # 2. Deterministic
  df_det = res_det %>%
    pivot_longer(cols = c(S, I, R), names_to = "State", values_to = "Count") %>%
    mutate(State = factor(State, levels = c("S", "I", "R")))
  
  # 3. plot
  p = ggplot() +

    geom_line(data = df_stoch, 
              aes(x = time, y = Count, group = interaction(sim_id, State), color = State), 
              alpha = alpha, linewidth = linewidth_stoch) +
    
    # scale_color_manual(values = c("S" = "steelblue", "I" = "firebrick", "R" = "chocolate1"),
    #                    guide = guide_legend(title = "Stochastic")) +
    
    scale_color_manual(values = c("S" = "blue", "I" = "red", "R" = "darkorange1"),
                       guide = guide_legend(title = "Status")) +
    
    
    # # Susceptible (steelblue)
    # geom_line(data = df_det %>% filter(State == "S"), 
    #           aes(x = time, y = Count), 
    #           color = "blue", alpha = 1, linewidth = linewidth_det) +
    # # Infected (firebrick)
    # geom_line(data = df_det %>% filter(State == "I"), 
    #           aes(x = time, y = Count), 
    #           color = "red", alpha = 1, linewidth = linewidth_det) +
    # # Recovered (chocolate1)
    # geom_line(data = df_det %>% filter(State == "R"), 
    #           aes(x = time, y = Count), 
    #           color = "darkorange1", alpha = 1, linewidth = linewidth_det) +
    
    
    
    geom_line(data = df_det %>% filter(State == "S"), 
              aes(x = time, y = Count), 
              color = "steelblue", alpha = 1, linewidth = linewidth_det) +
    # Infected (firebrick)
    geom_line(data = df_det %>% filter(State == "I"), 
              aes(x = time, y = Count), 
              color = "brown2", alpha = 1, linewidth = linewidth_det) +
    # Recovered (chocolate1)
    geom_line(data = df_det %>% filter(State == "R"), 
              aes(x = time, y = Count), 
              color = "chocolate1", alpha = 1, linewidth = linewidth_det) +
    
    
    scale_y_continuous(labels = scales::comma) + 
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time (days)",
      y = "Number of People"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}

# Deterministic SIR

sir_deterministic = function(N, I0, R0, d, period = 150) {
  
  # 1. Compute parameters
  gamma = 1 / d
  beta = R0 * gamma  # R0 = beta / gamma
  
  # 2. Initial State
  initial_state = c(
    S = N - I0, 
    I = I0, 
    R = 0
  )
  
  # 3. Differential Equations
  sir_equations = function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS = -beta * S * I / N
      dI =  beta * S * I / N - gamma * I
      dR =  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # 4. Setup
  params = c(beta = beta, gamma = gamma, N = N)
  times = seq(0, period, by = 1)
  
  # 5. Solve the equations
  out = ode(y = initial_state, times = times, func = sir_equations, parms = params)
  
  # 6. Convert into data.frame
  out_df = as.data.frame(out)
  
  return(out_df)
}

run_mean_field_sir = function(N = 10000, num_init_infected = 5, period = 150, R0 = 2.5, d = 7) {
  
  
  S_count = N - num_init_infected
  I_count = num_init_infected
  R_count = 0
  
  results = data.frame(time = 0, S = S_count, I = I_count, R = R_count)
  beta = R0 / d
  
  for (day in 1:period) {
    if (I_count == 0) {
      
      results = rbind(results, data.frame(time = day, S = S_count, I = 0, R = R_count))
      next
    }
    
    # 1. S -> I 
    prob_infection = 1 - exp(-beta * I_count / N)
    new_infections = rbinom(1, size = S_count, prob = prob_infection)
    
    # 2. I -> R
    prob_recovery = 1 / d
    new_recoveries = rbinom(1, size = I_count, prob = prob_recovery)
    
    # 3. Update the count
    S_count = S_count - new_infections
    I_count = I_count + new_infections - new_recoveries
    R_count = R_count + new_recoveries
    
    results = rbind(results, data.frame(time = day, S = S_count, I = I_count, R = R_count))
  }
  
  return(results)
}


generate_random_networks_empirical = function(N = 10000, target_k = 13.2, num_sims = 100, seed = 100) {
  
  set.seed(seed)
  er_adj_list = list()
  
  p_er = target_k / (N - 1)
  
  for (i in 1:num_sims) {
    g = igraph::sample_gnp(n = N, p = p_er, directed = FALSE)

    adj_matrix = igraph::as_adjacency_matrix(g, sparse = TRUE)
    
    er_adj_list[[i]] = adj_matrix
  }
  
  cat('Finished!\n')
  return(er_adj_list)
}

plot_sir_grid = function(sim_results_list, det_df, title = "SIR Model", R0 = 2.5, d = 7, period = 150, newpage = TRUE) {
  
  num_simulations = length(sim_results_list)
  N = det_df$S[1] + det_df$I[1] + det_df$R[1]
  
  gamma = 1 / d
  beta = R0 * gamma / N 
  
  if (newpage) grid.newpage()
  
    pushViewport(plotViewport(c(4, 4, 1, 1))) 
    pushViewport(dataViewport(xData = c(0, period), yData = c(0, N * 1.2), extension = 0))
    
      # Background / Axis / Scale
      grid.rect()
      grid.xaxis(at = seq(0, period, by = 30)) # 30 days
      grid.yaxis()
      
      # Title and Axis Label
      grid.text(paste0(title),
                y = unit(1, 'npc') - unit(0.75, 'lines'))
      grid.text('Time (days)', y = unit(-2.5, 'lines'))
      grid.text('Number of People', x = unit(-3.5, 'lines'), rot = 90)
      
      peak_stoch_vals = numeric(num_simulations)
      for (i in 1:num_simulations) {
        res = sim_results_list[[i]]
        grid.lines(x = res$time, y = res$S, default.units = "native", gp = gpar(col = 'blue', alpha = 0.1, lwd = 1.5))
        grid.lines(x = res$time, y = res$I, default.units = "native", gp = gpar(col = 'red', alpha = 0.1, lwd = 1.5))
        grid.lines(x = res$time, y = res$R, default.units = "native", gp = gpar(col = 'chocolate1', alpha = 0.1, lwd = 1.5))
        peak_stoch_vals[i] = max(res$I)
      }
      average_peak = mean(peak_stoch_vals)
      
      # --- Deterministic ---
      grid.lines(x = det_df$time, y = det_df$S, default.units = "native", gp = gpar(col = "darkblue", lwd = 2))
      grid.lines(x = det_df$time, y = det_df$I, default.units = "native", gp = gpar(col = "darkred", lwd = 2))
      grid.lines(x = det_df$time, y = det_df$R, default.units = "native", gp = gpar(col = "darkorange3", lwd = 2))
      
      theoretical_peak = max(det_df$I)
      peak_time_det = det_df$time[which.max(det_df$I)]
      
      # --- Peak info and dotted line ---
      grid.lines(x = c(0, period), y = c(theoretical_peak, theoretical_peak), default.units = "native", 
                 gp = gpar(lty = "dashed", col = "black", lwd = 1))
      
      peak_text = sprintf("Theoretical Peak: %d (Day %d)\nAvg. Stochastic Peak: %d", 
                          round(theoretical_peak), round(peak_time_det), round(average_peak))
      grid.text(peak_text, x = unit(0.95, "npc"), y = unit(theoretical_peak, "native") + unit(1, "lines"), 
                just = "right", gp = gpar(fontsize = 9))
      
      # --- Subtitle (params) ---
      param_expr = bquote(beta == .(sprintf("%.5f", beta)) ~ "|" ~ 
                          gamma == .(sprintf("%.4f", gamma)) ~ "|" ~ 
                          R[0] == .(sprintf("%.1f", R0)))
      
      grid.text(param_expr, 
                x = unit(0.5, 'npc'), 
                y = unit(1, 'npc') - unit(2, 'lines'), 
                gp = gpar(fontsize = 9))
      
    popViewport()
    
    # --- Legend ---
    pushViewport(viewport(x = 0.8, y = 0.7, width = 0.4, height = 0.2))
      grid.rect(gp = gpar(fill = 'white'))
      
      grid.text("Susceptible", x = 0.2, y = 0.8, just = "left", gp = gpar(fontsize = 9))
      grid.lines(x = c(0.05, 0.15), y = 0.8, gp = gpar(col = "blue", lwd = 2))
      
      grid.text("Infected", x = 0.2, y = 0.5, just = "left", gp = gpar(fontsize = 9))
      grid.lines(x = c(0.05, 0.15), y = 0.5, gp = gpar(col = "red", lwd = 2))
      
      grid.text("Recovered", x = 0.2, y = 0.2, just = "left", gp = gpar(fontsize = 9))
      grid.lines(x = c(0.05, 0.15), y = 0.2, gp = gpar(col = "orange", lwd = 2))
      
    popViewport() 
    popViewport() # dataViewport, plotViewport 
}

save_sir_grid_pdf = function(sim_results_list, det_df, filename, 
                             title = "SIR Model", R0 = 2.5, d = 7, period = 150,
                             width = 8, height = 8, pointsize = 10) {
  
  cm = 1 / 2.54
  pdf(file = filename, width = width * cm, height = height * cm, pointsize = pointsize)
  
  plot_sir_grid(
    sim_results_list = sim_results_list,
    det_df = det_df,
    title = title,
    R0 = R0,
    d = d,
    period = period
  )
  
  dev.off()
  
  cat(sprintf('PDF saved successfully: %s\n', filename))
}




save_sir_grid_combined_2x1 = function(
    res_1, det_df_1, title_1,  
    res_2, det_df_2, title_2,  
    filename, main_title = "Overall Title",
    R0 = 2.5, d = 7, period = 150,
    width = 16, height = 9.5, pointsize = 10
) {
  
  cm = 1 / 2.54
  pdf(file = filename, width = width * cm, height = height * cm, pointsize = pointsize)
  
  grid.newpage()

  lo = grid.layout(nrow = 2, ncol = 2, 
                   heights = unit(c(1.5, 1), c("cm", "null")),
                   widths = unit(c(1, 1), c("null", "null")))
  
  pushViewport(viewport(layout = lo))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(main_title, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  plot_sir_grid(res_1, det_df_1, title = title_1, R0 = R0, d = d, period = period, newpage = FALSE)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  plot_sir_grid(res_2, det_df_2, title = title_2, R0 = R0, d = d, period = period, newpage = FALSE)
  popViewport()
  
  popViewport()
  dev.off()
  
  cat(sprintf('Combined 2x1 PDF saved: %s\n', filename))
}

save_sir_grid_combined_2x2 = function(
    res_list, det_df_list, title_list, 
    filename, main_title = "Overall Title",
    R0 = 2.5, d = 7, period = 150,
    width = 16, height = 17.5, pointsize = 10
) {
  
  num_plots = length(res_list)
  
  cm = 1 / 2.54
  pdf(file = filename, width = width * cm, height = height * cm, pointsize = pointsize)
  
  grid.newpage()
  
  lo = grid.layout(nrow = 3, ncol = 2, 
                   heights = unit(c(1.5, 1, 1), c("cm", "null", "null")),
                   widths = unit(c(1, 1), c("null", "null")))
  
  pushViewport(viewport(layout = lo))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text(main_title, gp = gpar(fontsize = 14, fontface = "bold"))
  popViewport()
  
  for (i in 1:num_plots) {

    row_idx = ifelse(i <= 2, 2, 3) 
    col_idx = ifelse(i %% 2 != 0, 1, 2) 
    
    pushViewport(viewport(layout.pos.row = row_idx, layout.pos.col = col_idx))
    plot_sir_grid(res_list[[i]], det_df_list[[i]], title = title_list[[i]], 
                  R0 = R0, d = d, period = period, newpage = FALSE)
    popViewport()
  }
  
  popViewport()
  dev.off()
  
  cat(sprintf('Combined 2x2 PDF saved: %s\n', filename))
}



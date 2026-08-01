library(tidyverse)
library(grid)
library(Matrix)

create_sir_grob = function(sim_results_list, det_df, title = "SIR Model", R0 = 2.5, d = 7, period = 150) {
  
  num_simulations = length(sim_results_list)
  N = det_df$S[1] + det_df$I[1] + det_df$R[1]
  gamma = 1 / d
  beta = R0 * gamma / N 
  
  main_vp = plotViewport(margins = c(4, 4, 1, 1), name = "main_vp")
  
  title_grob = textGrob(title, y = unit(1, 'npc') - unit(0.75, 'lines'))
  xlab_grob = textGrob('Time (days)', y = unit(-2.5, 'lines'))

  ylab_grob = textGrob('Number of People (K)', x = unit(-3, 'lines'), rot = 90)
  

  data_vp = dataViewport(xData = c(0, period), yData = c(0, N * 1.2), extension = 0, name = "data_vp")
  
  bg = rectGrob()
  x_axis = xaxisGrob(at = seq(0, period, by = 30))
  
  y_breaks = pretty(c(0, N * 1.2))      
  y_labels = y_breaks / 1000            
  y_axis = yaxisGrob(at = y_breaks, label = y_labels)
  
  stoch_grobs = list()
  peak_stoch_vals = numeric(num_simulations)
  for (i in 1:num_simulations) {
    res = sim_results_list[[i]]
    stoch_grobs[[paste0("S_", i)]] = linesGrob(x = res$time, y = res$S, default.units = "native", gp = gpar(col = 'blue', alpha = 0.1, lwd = 1.5))
    stoch_grobs[[paste0("I_", i)]] = linesGrob(x = res$time, y = res$I, default.units = "native", gp = gpar(col = 'red', alpha = 0.1, lwd = 1.5))
    stoch_grobs[[paste0("R_", i)]] = linesGrob(x = res$time, y = res$R, default.units = "native", gp = gpar(col = 'chocolate1', alpha = 0.1, lwd = 1.5))
    peak_stoch_vals[i] = max(res$I)
  }
  average_peak = mean(peak_stoch_vals)
  stoch_tree = do.call(gList, stoch_grobs)
  
  det_S = linesGrob(x = det_df$time, y = det_df$S, default.units = "native", gp = gpar(col = "darkblue", lwd = 2))
  det_I = linesGrob(x = det_df$time, y = det_df$I, default.units = "native", gp = gpar(col = "darkred", lwd = 2))
  det_R = linesGrob(x = det_df$time, y = det_df$R, default.units = "native", gp = gpar(col = "darkorange3", lwd = 2))
  
  deterministic_peak = max(det_df$I)
  peak_line = linesGrob(x = c(0, period), y = c(deterministic_peak, deterministic_peak), default.units = "native", gp = gpar(lty = "dashed", col = "black", lwd = 1))
  
  data_tree = gTree(
    children = gList(bg, x_axis, y_axis, stoch_tree, det_S, det_I, det_R, peak_line),
    vp = data_vp
  )
  
  final_grob = gTree(
    children = gList(title_grob, xlab_grob, ylab_grob, data_tree),
    vp = main_vp,
    name = paste0("sir_grob_", gsub(" ", "_", title))
  )
  
  return(final_grob)
}



create_sir_legend_grob = function() {
  
  leg_bg = rectGrob(gp = gpar(fill = 'white', col = 'black', lwd = 1))
  
  leg_t1 = textGrob("Susceptible", x = 0.35, y = 0.8, just = "left", gp = gpar(fontsize = 9))
  leg_l1 = linesGrob(x = c(0.1, 0.25), y = 0.8, gp = gpar(col = "blue", lwd = 2))
  
  leg_t2 = textGrob("Infected", x = 0.35, y = 0.5, just = "left", gp = gpar(fontsize = 9))
  leg_l2 = linesGrob(x = c(0.1, 0.25), y = 0.5, gp = gpar(col = "red", lwd = 2))
  
  leg_t3 = textGrob("Recovered", x = 0.35, y = 0.2, just = "left", gp = gpar(fontsize = 9))
  leg_l3 = linesGrob(x = c(0.1, 0.25), y = 0.2, gp = gpar(col = "chocolate1", lwd = 2))
  
  legend_tree = gTree(
    children = gList(leg_bg, leg_t1, leg_l1, leg_t2, leg_l2, leg_t3, leg_l3),
    name = "sir_legend_grob"
  )
  
  return(legend_tree)
}
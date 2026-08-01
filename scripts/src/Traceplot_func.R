library(grid)
library(tidyverse)

prep_multi_chain = function(mcmc_results_list) {
  
  num_stages = length(mcmc_results_list)
  
  df_list = lapply(1:num_stages, function(i) {
    res = mcmc_results_list[[i]]
    
    chain_df = as.data.frame(res$chain)
    
    inf_val = if (!is.null(res$config$inflation_scale)) {
      res$config$inflation_scale
    } else {
      NA
    }
    
    chain_df %>%
      mutate(
        stage = i,
        inflation = inf_val
      )
  })
  
  base_df = bind_rows(df_list) %>%
    mutate(global_iter = row_number())
  
  connections = base_df %>%
    group_by(stage) %>%
    slice_tail(n = 1) %>%        
    ungroup() %>%
    filter(stage < num_stages) %>% 
    mutate(stage = stage + 1)      
  
  final_df = bind_rows(base_df, connections) %>%
    arrange(global_iter, stage)
  
  return(final_df)
}


prep_single_chain = function(mcmc_result) {
  
  chain_df = as.data.frame(mcmc_result$chain)
  
  inf_val = if (!is.null(mcmc_result$config$inflation_scale)) {
    mcmc_result$config$inflation_scale
  } else {
    NA
  }
  
  final_df = chain_df %>%
    mutate(
      stage = 1,
      inflation = inf_val,
      global_iter = row_number()
    )
  
  return(final_df)
}


create_trace_grob = function(data, param_col, param_label = param_col, colors = NULL, main_title = NULL) {
  
  x_val = data$global_iter
  y_val = data[[param_col]]
  x_scale = c(min(x_val), max(x_val))
  y_range = max(y_val) - min(y_val)
  if(y_range == 0) y_range = 1e-6 
  y_scale = c(min(y_val) - y_range * 0.05, max(y_val) + y_range * 0.25)
  
  stages = unique(data$stage)
  num_stages = length(stages)
  if (is.null(colors)) {
    if (num_stages == 1) {
      colors = "steelblue"
    } else {
      custom_palette = colorRampPalette(c("navy", "steelblue", "darkorchid", "violetred", "firebrick"))
      colors = custom_palette(num_stages)
    }
  }
  
  vp = plotViewport(margins = c(4, 4, 1, 1), 
                    xscale = x_scale, yscale = y_scale,
                    name = paste0("vp_", param_col))
  
  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  
  line_grobs = list()
  for (s in stages) {
    stage_data = data[data$stage == s, ]
    line_grobs[[s]] = linesGrob(
      x = stage_data$global_iter, y = stage_data[[param_col]], 
      default.units = "native", gp = gpar(col = colors[s], lwd = 0.8) 
    )
  }
  lines_tree = do.call(gList, line_grobs)
  
  bound_grobs = list()
  if (num_stages > 1) {
    for (s in 1:(num_stages - 1)) {
      bound_x = max(data$global_iter[data$stage == s])
      bound_grobs[[s]] = linesGrob(
        x = unit(c(bound_x, bound_x), "native"),
        y = unit(c(0, 0.85), "npc"), 
        gp = gpar(col = "gray40", lty = "dashed", lwd = 1)
      )
    }
  }
  bounds_tree = do.call(gList, bound_grobs)
  
  x_breaks = pretty(x_scale)
  x_labels = x_breaks / 1000
  x_axis = xaxisGrob(at = x_breaks, label = x_labels)
  y_axis = yaxisGrob()
  x_lab = textGrob("Iteration (K)", y = unit(-2.5, "lines"))
  
  display_label = if (!is.null(main_title)) main_title else param_label
  
  label_cex = if (!is.null(main_title)) 1.0 else 1.25
  
  inner_param_lab = textGrob(
    display_label, 
    x = unit(0.5, "npc"),                    
    y = unit(1, "npc") - unit(1, "lines"),    
    just = c("center", "center"), 
    gp = gpar(cex = label_cex)
  )
  
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  result_grob = gTree(
    children = gList(bg, bounds_tree, lines_tree, x_axis, y_axis, x_lab, inner_param_lab, border_box),
    vp = vp,
    name = paste0("trace_", param_col)
  )
  
  return(result_grob)
}


param_dict_M2 = list(
  alpha          = expression(alpha),
  gamma_school   = expression(gamma[school]),
  gamma_pri_low  = expression(gamma[pri_low]),
  gamma_pri_high = expression(gamma[pri_high]),
  gamma_sec      = expression(gamma[sec]),
  gamma_ter      = expression(gamma[ter]),
  w              = expression(w),
  lambda         = expression(lambda),
  sigma_hhd      = expression(sigma[hhd]),
  mu             = expression(mu),
  phi_mm         = expression(phi[mm]),
  phi_ff         = expression(phi[ff]),
  sigma_re       = expression(sigma[re]),
  theta          = expression(theta),
  sigma          = expression(sigma)
)


generate_trace_grobs = function(data, target_params, dict = param_dict_M2, colors = NULL, main_title = NULL) {
  
  grobs_list = list()
  
  for (param in target_params) {
    p_label = if (!is.null(dict[[param]])) dict[[param]] else param
    
    grobs_list[[param]] = create_trace_grob(
      data = data, 
      param_col = param, 
      param_label = p_label, 
      colors = colors,
      main_title = main_title 
    )
  }
  
  return(grobs_list)
}





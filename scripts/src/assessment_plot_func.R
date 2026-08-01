library(grid)
library(tidyverse)

create_mean_degree_grob_simple = function(data, model_name, is_polymod_panel = FALSE) {
  
  poly_summary = data %>% filter(source == "POLYMOD")
  sim_summary  = data %>% filter(source == "Simulation")
  
  if (is_polymod_panel) {
    plot_data_main = poly_summary
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = sim_summary
    plot_data_bg   = poly_summary
    color_main     = "steelblue"
    color_bg       = "palevioletred1" 
  }
  
  old_labels = levels(poly_summary$age_group)
  num_groups = length(old_labels)
  
  breaks_5yr = seq(0, num_groups * 5, by = 5)
  new_labels = paste0("[", breaks_5yr[-length(breaks_5yr)], ", ", breaks_5yr[-1], ")")
  new_labels[num_groups] = paste0("[", breaks_5yr[num_groups], ", ", breaks_5yr[num_groups + 1], "]")
  
  x_scale = c(1, num_groups)
  y_scale = c(0, 24) 
  
  main_vp = plotViewport(margins = c(5, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")
  
  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  x_axis_line = linesGrob(x = unit(c(1, num_groups), "native"), y = unit(0, "npc"), gp = gpar(col = "black"))
  x_ticks = segmentsGrob(
    x0 = unit(1:num_groups, "native"), y0 = unit(0, "npc"),
    x1 = unit(1:num_groups, "native"), y1 = unit(-1.5, "mm"),
    gp = gpar(col = "black")
  )
  x_axis = gTree(children = gList(x_axis_line, x_ticks))
  
  x_labels_grob = textGrob(
    label = new_labels, 
    x = unit(1:num_groups, "native"),
    y = unit(-2.5, "mm"),
    just = c("right", "center"), 
    rot = 90, 
    gp = gpar(fontsize = 8)
  )
  
  y_axis = yaxisGrob()
  
  xlab_grob = textGrob("Age Group", y = unit(-3.8, "lines")) 
  ylab_grob = textGrob("Mean Total Contacts", x = unit(-3, "lines"), rot = 90)
  
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  
  lines_list = list()
  
  if (!is.null(plot_data_bg)) {
    bg_lines = linesGrob(x = as.numeric(plot_data_bg$age_group), y = plot_data_bg$mean_degree, 
                         default.units = "native", gp = gpar(col = color_bg, lwd = 2, lty = "solid"))
    bg_points = pointsGrob(x = as.numeric(plot_data_bg$age_group), y = plot_data_bg$mean_degree, 
                           pch = 16, size = unit(0.5, "char"), default.units = "native", gp = gpar(col = color_bg))
    lines_list = append(lines_list, list(bg_lines, bg_points))
  }
  
  main_lines = linesGrob(x = as.numeric(plot_data_main$age_group), y = plot_data_main$mean_degree, 
                         default.units = "native", gp = gpar(col = color_main, lwd = 2.5))
  main_points = pointsGrob(x = as.numeric(plot_data_main$age_group), y = plot_data_main$mean_degree, 
                           pch = 16, size = unit(0.5, "char"), default.units = "native", gp = gpar(col = color_main))
  lines_list = append(lines_list, list(main_lines, main_points))
  
  data_tree = gTree(
    children = gList(bg, x_axis, x_labels_grob, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  final_grob = gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp)
  
  return(final_grob)
}

create_shared_legend_grob = function() {
  
  leg_bg = rectGrob(gp = gpar(fill = "white", col = NA))
  
  leg_poly_line = linesGrob(x = c(0.1, 0.3), y = 0.7, gp = gpar(col = "red", lwd = 2.5))
  leg_poly_text = textGrob("POLYMOD", x = 0.35, y = 0.7, just = "left", gp = gpar(fontsize = 10))
  
  leg_sim_line = linesGrob(x = c(0.1, 0.3), y = 0.3, gp = gpar(col = "steelblue", lwd = 2.5))
  leg_sim_text = textGrob("Model", x = 0.35, y = 0.3, just = "left", gp = gpar(fontsize = 10))
  
  legend_tree = gTree(
    children = gList(leg_bg, leg_poly_line, leg_poly_text, leg_sim_line, leg_sim_text)
  )
  
  return(legend_tree)
}

create_degree_cdf_grob_simple = function(data, model_name, is_polymod_panel = FALSE, max_degree = 40) {
  
  data = data %>% filter(degree <= max_degree)
  
  prepare_step_data = function(df) {
    df = df %>% arrange(degree)
    total = sum(df$mean_count)
    df = df %>% mutate(cum_prob = cumsum(mean_count) / total)
    
    n = nrow(df)
    step_x = numeric(n * 2)
    step_y = numeric(n * 2)
    
    step_x[1] = df$degree[1]
    step_y[1] = df$cum_prob[1]
    
    for (i in 1:(n - 1)) {
      step_x[2 * i] = df$degree[i + 1]
      step_y[2 * i] = df$cum_prob[i]
      step_x[2 * i + 1] = df$degree[i + 1]
      step_y[2 * i + 1] = df$cum_prob[i + 1]
    }
    step_x[n * 2] = max_degree
    step_y[n * 2] = df$cum_prob[n]
    
    return(data.frame(x = step_x, y = step_y))
  }
  
  poly_summary = data %>% filter(source == "POLYMOD")
  sim_summary  = data %>% filter(source == "Simulation")
  
  if (is_polymod_panel) {
    plot_data_main = prepare_step_data(poly_summary)
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = prepare_step_data(sim_summary)
    plot_data_bg   = prepare_step_data(poly_summary)
    color_main     = "steelblue"
    color_bg       = "palevioletred1"
  }
  
  x_scale = c(0, max_degree)
  y_scale = c(0, 1.2)
  
  main_vp = plotViewport(margins = c(4, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")

  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  x_ticks_at = seq(0, max_degree, by = 5)
  x_labels_str = ifelse(x_ticks_at %% 10 == 0, as.character(x_ticks_at), "")
  
  x_axis = xaxisGrob(at = x_ticks_at, label = x_labels_str, gp = gpar(fontsize = 10))
  
  y_axis = yaxisGrob(at = seq(0, 1, by = 0.2))
  
  xlab_grob = textGrob("Number of Contacts (Degree)", y = unit(-2.5, "lines"))
  ylab_grob = textGrob("Cumulative Probability", x = unit(-3, "lines"), rot = 90)
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  
  lines_list = list()
  if (!is.null(plot_data_bg)) {
    lines_list = append(lines_list, list(linesGrob(x = plot_data_bg$x, y = plot_data_bg$y, 
                                                   default.units = "native", gp = gpar(col = color_bg, lwd = 2))))
  }
  lines_list = append(lines_list, list(linesGrob(x = plot_data_main$x, y = plot_data_main$y, 
                                                 default.units = "native", gp = gpar(col = color_main, lwd = 2))))
  
  data_tree = gTree(
    children = gList(bg, x_axis, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  return(gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp))
}

create_agediff_cdf_grob_simple = function(data_sim, data_poly, model_name, is_polymod_panel = FALSE, max_diff = 80) {
  
  prepare_step_data = function(df) {
    df_counts = df %>%
      filter(abs_age_diff <= max_diff) %>%
      count(abs_age_diff, name = "count") %>%
      arrange(abs_age_diff)
    
    total = sum(df_counts$count)
    df_counts = df_counts %>% mutate(cum_prob = cumsum(count) / total)
    
    n = nrow(df_counts)
    step_x = numeric(n * 2)
    step_y = numeric(n * 2)
    
    step_x[1] = df_counts$abs_age_diff[1]
    step_y[1] = df_counts$cum_prob[1]
    
    if (n > 1) {
      for (i in 1:(n - 1)) {
        step_x[2 * i] = df_counts$abs_age_diff[i + 1]
        step_y[2 * i] = df_counts$cum_prob[i]
        step_x[2 * i + 1] = df_counts$abs_age_diff[i + 1]
        step_y[2 * i + 1] = df_counts$cum_prob[i + 1]
      }
    }
    step_x[n * 2] = max_diff
    step_y[n * 2] = df_counts$cum_prob[n]
    
    return(data.frame(x = step_x, y = step_y))
  }
  
  if (is_polymod_panel) {
    plot_data_main = prepare_step_data(data_poly)
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = prepare_step_data(data_sim)
    plot_data_bg   = prepare_step_data(data_poly)
    color_main     = "steelblue"
    color_bg       = "palevioletred1"
  }
  
  x_scale = c(0, max_diff)
  y_scale = c(0, 1.2) 
  
  main_vp = plotViewport(margins = c(4, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")
  
  
  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  x_ticks_at = seq(0, max_diff, by = 5)
  x_labels_str = ifelse(x_ticks_at %% 10 == 0, as.character(x_ticks_at), "")
  
  x_axis = xaxisGrob(at = x_ticks_at, label = x_labels_str, gp = gpar(fontsize = 10))
  
  y_axis = yaxisGrob(at = seq(0, 1, by = 0.2))
  
  xlab_grob = textGrob("Absolute Age Difference (Years)", y = unit(-2.5, "lines"))
  ylab_grob = textGrob("Cumulative Probability", x = unit(-3, "lines"), rot = 90)
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  
  lines_list = list()
  if (!is.null(plot_data_bg)) {
    lines_list = append(lines_list, list(linesGrob(x = plot_data_bg$x, y = plot_data_bg$y, 
                                                   default.units = "native", gp = gpar(col = color_bg, lwd = 2))))
  }
  lines_list = append(lines_list, list(linesGrob(x = plot_data_main$x, y = plot_data_main$y, 
                                                 default.units = "native", gp = gpar(col = color_main, lwd = 2))))
  
  data_tree = gTree(
    children = gList(bg, x_axis, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  return(gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp))
}


create_diag_1yr_grob_simple = function(data, model_name, is_polymod_panel = FALSE, max_age = 80, max_y = 6) {
  
  poly_summary = data %>% filter(source == "POLYMOD", age <= max_age)
  sim_summary  = data %>% filter(source == "Simulation", age <= max_age)
  
  if (is_polymod_panel) {
    plot_data_main = poly_summary
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = sim_summary
    plot_data_bg   = poly_summary
    color_main     = "steelblue"
    color_bg       = "palevioletred1" 
  }
  
  x_scale = c(0, max_age)
  y_scale = c(0, max_y)
  
  main_vp = plotViewport(margins = c(4, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")
  
  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))

  x_ticks_at = seq(0, max_age, by = 5)
  x_labels_str = ifelse(x_ticks_at %% 10 == 0, as.character(x_ticks_at), "")
  x_axis = xaxisGrob(at = x_ticks_at, label = x_labels_str, gp = gpar(fontsize = 10))
  
  y_axis = yaxisGrob()
  
  xlab_grob = textGrob("Age (Years)", y = unit(-2.5, "lines"))
  ylab_grob = textGrob("Mean Diagonal Contacts", x = unit(-3, "lines"), rot = 90)
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  
  lines_list = list()
  
  if (!is.null(plot_data_bg)) {
    bg_lines = linesGrob(x = plot_data_bg$age, y = plot_data_bg$mean_contacts, 
                         default.units = "native", gp = gpar(col = color_bg, lwd = 1.5))
    lines_list = append(lines_list, list(bg_lines))
  }
  

  main_lines = linesGrob(x = plot_data_main$age, y = plot_data_main$mean_contacts, 
                         default.units = "native", gp = gpar(col = color_main, lwd = 1.5))
  lines_list = append(lines_list, list(main_lines))
  
  data_tree = gTree(
    children = gList(bg, x_axis, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  final_grob = gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp)
  
  return(final_grob)
}


create_agediff_dist_grob_simple = function(data_sim, data_poly, model_name, is_polymod_panel = FALSE, max_diff = 80, max_y = 0.15) {
  
  prepare_dist_data = function(df) {
    df_counts = df %>%
      filter(abs_age_diff <= max_diff) %>%
      count(abs_age_diff, name = "count") %>%
      arrange(abs_age_diff)
    
    total = sum(df_counts$count)
    df_counts = df_counts %>% mutate(prob = count / total)
    
    return(df_counts)
  }
  
  if (is_polymod_panel) {
    plot_data_main = prepare_dist_data(data_poly)
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = prepare_dist_data(data_sim)
    plot_data_bg   = prepare_dist_data(data_poly)
    color_main     = "steelblue"
    color_bg       = "palevioletred1"
  }
  

  x_scale = c(0, max_diff)
  y_scale = c(0, max_y) 
  
  main_vp = plotViewport(margins = c(4, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")
  

  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  

  x_ticks_at = seq(0, max_diff, by = 5)
  x_labels_str = ifelse(x_ticks_at %% 10 == 0, as.character(x_ticks_at), "")
  x_axis = xaxisGrob(at = x_ticks_at, label = x_labels_str, gp = gpar(fontsize = 10))
  
  y_axis = yaxisGrob()
  
  xlab_grob = textGrob("Absolute Age Difference (Years)", y = unit(-2.5, "lines"))
  ylab_grob = textGrob("Probability", x = unit(-3.5, "lines"), rot = 90)
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  

  lines_list = list()
  if (!is.null(plot_data_bg)) {
    lines_list = append(lines_list, list(linesGrob(x = plot_data_bg$abs_age_diff, y = plot_data_bg$prob, 
                                                   default.units = "native", gp = gpar(col = color_bg, lwd = 2))))
  }
  lines_list = append(lines_list, list(linesGrob(x = plot_data_main$abs_age_diff, y = plot_data_main$prob, 
                                                 default.units = "native", gp = gpar(col = color_main, lwd = 2.5))))
  
  data_tree = gTree(
    children = gList(bg, x_axis, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  return(gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp))
}


create_diag_step_grob_simple = function(data, model_name, is_polymod_panel = FALSE, max_y = 10) {
  
  poly_summary = data %>% filter(source == "POLYMOD")
  sim_summary  = data %>% filter(source == "Simulation")
  
  if (is_polymod_panel) {
    plot_data_main = poly_summary
    plot_data_bg   = NULL
    color_main     = "red"
  } else {
    plot_data_main = sim_summary
    plot_data_bg   = poly_summary
    color_main     = "steelblue"
    color_bg       = "palevioletred1" 
  }
  
  age_groups = levels(poly_summary$age_group)
  num_groups = length(age_groups)
  
  breaks_5yr = seq(0, num_groups * 5, by = 5)
  new_labels = paste0("[", breaks_5yr[-length(breaks_5yr)], ", ", breaks_5yr[-1], ")")
  new_labels[num_groups] = paste0("[", breaks_5yr[num_groups], ", ", breaks_5yr[num_groups + 1], "]")
  
  x_scale = c(1, num_groups)
  y_scale = c(0, max_y)
  
  main_vp = plotViewport(margins = c(5, 4, 1, 1), name = "main_vp")
  data_vp = dataViewport(xData = x_scale, yData = y_scale, extension = 0, name = "data_vp")
  

  bg = rectGrob(gp = gpar(fill = "white", col = NA))
  border_box = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  

  x_axis_line = linesGrob(x = unit(c(1, num_groups), "native"), y = unit(0, "npc"))
  x_ticks = segmentsGrob(
    x0 = unit(1:num_groups, "native"), y0 = unit(0, "npc"),
    x1 = unit(1:num_groups, "native"), y1 = unit(-1.5, "mm")
  )
  x_axis = gTree(children = gList(x_axis_line, x_ticks))
  

  x_labels_grob = textGrob(
    label = new_labels, 
    x = unit(1:num_groups, "native"), y = unit(-2.5, "mm"), 
    just = c("right", "center"), rot = 90, gp = gpar(fontsize = 8)
  )
  
  y_axis = yaxisGrob()
  
  xlab_grob = textGrob("Age Group", y = unit(-3.8, "lines")) 
  ylab_grob = textGrob("Mean Diagonal Contacts", x = unit(-3, "lines"), rot = 90)
  title_grob = textGrob(model_name, x = 0.5, y = unit(1, "npc") - unit(1, "lines"), just = "center")
  
  lines_list = list()
  if (!is.null(plot_data_bg)) {
    bg_lines = linesGrob(x = as.numeric(plot_data_bg$age_group), y = plot_data_bg$mean_contacts, 
                         default.units = "native", gp = gpar(col = color_bg, lwd = 2))
    bg_points = pointsGrob(x = as.numeric(plot_data_bg$age_group), y = plot_data_bg$mean_contacts, 
                           pch = 16, size = unit(0.5, "char"), default.units = "native", gp = gpar(col = color_bg))
    lines_list = append(lines_list, list(bg_lines, bg_points))
  }
  
  main_lines = linesGrob(x = as.numeric(plot_data_main$age_group), y = plot_data_main$mean_contacts, 
                         default.units = "native", gp = gpar(col = color_main, lwd = 2.5))
  main_points = pointsGrob(x = as.numeric(plot_data_main$age_group), y = plot_data_main$mean_contacts, 
                           pch = 16, size = unit(0.5, "char"), default.units = "native", gp = gpar(col = color_main))
  lines_list = append(lines_list, list(main_lines, main_points))
  
  data_tree = gTree(
    children = gList(bg, x_axis, x_labels_grob, y_axis, do.call(gList, lines_list), border_box, title_grob),
    vp = data_vp
  )
  return(gTree(children = gList(xlab_grob, ylab_grob, data_tree), vp = main_vp))
}


create_heatmap_grob_simple = function(data, model_name, is_polymod_panel = FALSE, max_val = 5.0) {
  
  target_data = data %>% filter(source == (if(is_polymod_panel) "POLYMOD" else "Simulation"))
  

  breaks_5yr = seq(0, 80, by = 5)
  labels_5yr = paste0(breaks_5yr[-length(breaks_5yr)], "-", breaks_5yr[-1] - 1)
  labels_5yr[length(labels_5yr)] = "75-80"
  
  target_data = target_data %>%
    mutate(
      ego_group = factor(ego_group, levels = labels_5yr),
      alter_group = factor(alter_group, levels = labels_5yr)
    )
  
  age_groups = levels(target_data$ego_group)
  n = length(age_groups) 
  
  blue_pal = colorRampPalette(c("white", "skyblue2", "dodgerblue3", "royalblue3","blue3","navy"))(100)

  target_data = target_data %>%
    mutate(
      intensity = pmin(mean_contacts / max_val, 1),
      col_idx = as.integer(intensity * 99) + 1,
      fill_col = blue_pal[col_idx]
    )


  main_vp = plotViewport(margins = c(3, 4, 2, 1), name = "main_vp")
  data_vp = viewport(xscale = c(0.5, n + 0.5), yscale = c(0.5, n + 0.5), name = "data_vp")
  
  cells = rectGrob(
    x = unit(as.numeric(target_data$alter_group), "native"),
    y = unit(as.numeric(target_data$ego_group), "native"),
    width = unit(1, "native"), height = unit(1, "native"),
    gp = gpar(fill = target_data$fill_col, col = NA)
  )
  
  border = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  major_ticks_at = seq(0.5, n + 0.5, by = 2) # 0, 10, 20... 80 
  minor_ticks_at = seq(1.5, n - 0.5, by = 2) # 5, 15, 25... 75 
  major_labels_str = as.character(seq(0, 80, by = 10))
  
  x_ticks_major = segmentsGrob(
    x0 = unit(major_ticks_at, "native"), y0 = unit(0, "npc"),
    x1 = unit(major_ticks_at, "native"), y1 = unit(-2.0, "mm"), gp = gpar(col = "black")
  )
  x_ticks_minor = segmentsGrob(
    x0 = unit(minor_ticks_at, "native"), y0 = unit(0, "npc"),
    x1 = unit(minor_ticks_at, "native"), y1 = unit(-1.5, "mm"), gp = gpar(col = "black")
  )
  x_labels_grob = textGrob(
    label = major_labels_str, 
    x = unit(major_ticks_at, "native"), 
    y = unit(-3.0, "mm"), just = "top", gp = gpar(fontsize = 10)
  )
  x_axis = gTree(children = gList(x_ticks_major, x_ticks_minor, x_labels_grob))
  
  y_ticks_major = segmentsGrob(
    y0 = unit(major_ticks_at, "native"), x0 = unit(0, "npc"),
    y1 = unit(major_ticks_at, "native"), x1 = unit(-2.0, "mm"), gp = gpar(col = "black")
  )
  y_ticks_minor = segmentsGrob(
    y0 = unit(minor_ticks_at, "native"), x0 = unit(0, "npc"),
    y1 = unit(minor_ticks_at, "native"), x1 = unit(-1.5, "mm"), gp = gpar(col = "black")
  )
  y_labels_grob = textGrob(
    label = major_labels_str, 
    y = unit(major_ticks_at, "native"), 
    x = unit(-3.0, "mm"), just = "right", gp = gpar(fontsize = 10)
  )
  y_axis = gTree(children = gList(y_ticks_major, y_ticks_minor, y_labels_grob))
  # -----------------------------------
  
  xlab = textGrob("Contact's Age (Years)", y = unit(-2.25, "lines"))
  ylab = textGrob("Participant's Age (Years)", x = unit(-3, "lines"), rot = 90)
  title = textGrob(model_name, y = unit(1, "npc") + unit(0.6, "lines"))
  
  data_tree = gTree(children = gList(cells, border, x_axis, y_axis, title), vp = data_vp)
  return(gTree(children = gList(xlab, ylab, data_tree), vp = main_vp))
}


create_residual_heatmap_grob_simple = function(data, model_name, max_res = 2.0) {
  
  poly_data = data %>% filter(source == "POLYMOD") %>% select(ego_group, alter_group, poly_contacts = mean_contacts)
  sim_data  = data %>% filter(source == "Simulation") %>% select(ego_group, alter_group, sim_contacts = mean_contacts)
  
  breaks_5yr = seq(0, 80, by = 5)
  labels_5yr = paste0(breaks_5yr[-length(breaks_5yr)], "-", breaks_5yr[-1] - 1)
  labels_5yr[length(labels_5yr)] = "75-80"
  
  target_data = poly_data %>%
    full_join(sim_data, by = c("ego_group", "alter_group")) %>%
    mutate(
      poly_contacts = replace_na(poly_contacts, 0), 
      sim_contacts = replace_na(sim_contacts, 0),   
      ego_group = factor(ego_group, levels = labels_5yr),
      alter_group = factor(alter_group, levels = labels_5yr),
      residual = poly_contacts - sim_contacts,
      norm_res = pmax(pmin(residual / max_res, 1), -1),
      idx = (norm_res + 1) / 2
    ) %>%
    filter(!is.na(idx))
  
  n = length(levels(target_data$ego_group))
  
  div_pal = colorRamp(c("#0571B0", "white", "#CA0020"))
  rgb_matrix = div_pal(target_data$idx) / 255
  target_data$fill_col = rgb(rgb_matrix[,1], rgb_matrix[,2], rgb_matrix[,3])
  
  main_vp = plotViewport(margins = c(3, 4, 2, 1), name = "main_vp")
  data_vp = viewport(xscale = c(0.5, n + 0.5), yscale = c(0.5, n + 0.5), name = "data_vp")
  
  cells = rectGrob(
    x = unit(as.numeric(target_data$alter_group), "native"),
    y = unit(as.numeric(target_data$ego_group), "native"),
    width = unit(1, "native"), height = unit(1, "native"),
    gp = gpar(fill = target_data$fill_col, col = NA)
  )
  
  border = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  major_ticks_at = seq(0.5, n + 0.5, by = 2) 
  minor_ticks_at = seq(1.5, n - 0.5, by = 2) 
  major_labels_str = as.character(seq(0, 80, by = 10))
  
  x_ticks_major = segmentsGrob(
    x0 = unit(major_ticks_at, "native"), y0 = unit(0, "npc"),
    x1 = unit(major_ticks_at, "native"), y1 = unit(-2.0, "mm"), gp = gpar(col = "black")
  )
  x_ticks_minor = segmentsGrob(
    x0 = unit(minor_ticks_at, "native"), y0 = unit(0, "npc"),
    x1 = unit(minor_ticks_at, "native"), y1 = unit(-1.5, "mm"), gp = gpar(col = "black")
  )
  x_labels_grob = textGrob(
    label = major_labels_str, 
    x = unit(major_ticks_at, "native"), 
    y = unit(-3.0, "mm"), just = "top", gp = gpar(fontsize = 10)
  )
  x_axis = gTree(children = gList(x_ticks_major, x_ticks_minor, x_labels_grob))
  
  y_ticks_major = segmentsGrob(
    y0 = unit(major_ticks_at, "native"), x0 = unit(0, "npc"),
    y1 = unit(major_ticks_at, "native"), x1 = unit(-2.0, "mm"), gp = gpar(col = "black")
  )
  y_ticks_minor = segmentsGrob(
    y0 = unit(minor_ticks_at, "native"), x0 = unit(0, "npc"),
    y1 = unit(minor_ticks_at, "native"), x1 = unit(-1.5, "mm"), gp = gpar(col = "black")
  )
  y_labels_grob = textGrob(
    label = major_labels_str, 
    y = unit(major_ticks_at, "native"), 
    x = unit(-3.0, "mm"), just = "right", gp = gpar(fontsize = 10)
  )
  y_axis = gTree(children = gList(y_ticks_major, y_ticks_minor, y_labels_grob))
  
  xlab = textGrob("Contact's Age (Years)", y = unit(-2.25, "lines"))
  ylab = textGrob("Participant's Age (Years)", x = unit(-3, "lines"), rot = 90)
  
  title_text = paste0(model_name, " Residual")
  title = textGrob(title_text, y = unit(1, "npc") + unit(0.6, "lines"))
  
  data_tree = gTree(children = gList(cells, border, x_axis, y_axis, title), vp = data_vp)
  return(gTree(children = gList(xlab, ylab, data_tree), vp = main_vp))
}



create_heatmap_legend_grob = function(max_val = 10) {
  # Define the same color palette as the main heatmap
  blue_pal = colorRampPalette(c("white", "skyblue2", "dodgerblue3", "royalblue3","blue3","navy"))(100)
  
  # Create a matrix to draw the color gradient
  grad_matrix = matrix(rev(blue_pal), ncol = 1) 
  
  # Reduce height to 0.4 (40%) and shift left (x=0.2) to prevent cropping
  # Fix the width of the color bar itself to 5mm
  legend_vp = viewport(x = 0.2, y = 0.5, width = unit(5, "mm"), height = 0.4, just = "left")
  
  # Color bar body
  color_bar = rasterGrob(grad_matrix, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  border = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  # [Modification] Setup ticks and labels with 4 intervals (e.g., 0, 2.5, 5, 7.5, 10)
  y_ticks_at = seq(0, 1, by = 0.25) 
  labels_str = as.character(seq(0, max_val, length.out = 5))
  
  # Optional: Format labels to have consistent decimal places if needed (e.g., "2.5", "5.0")
  # labels_str = sprintf("%.1f", seq(0, max_val, length.out = 5))
  
  ticks = segmentsGrob(
    x0 = unit(1, "npc"), y0 = unit(y_ticks_at, "npc"),
    x1 = unit(1, "npc") + unit(2, "mm"), y1 = unit(y_ticks_at, "npc"),
    gp = gpar(col = "black")
  )
  
  labels = textGrob(
    label = labels_str, 
    x = unit(1, "npc") + unit(3, "mm"), y = unit(y_ticks_at, "npc"), 
    just = "left", gp = gpar(fontsize = 10)
  )
  
  # Shift the title slightly upwards (+2.5 lines)
  title = textGrob("Mean\nContacts", 
                   x = unit(0.5, "npc"), 
                   y = unit(1, "npc") + unit(2.5, "lines"), 
                   just = "center", gp = gpar(fontsize = 10, lineheight = 0.9))
  
  # Combine all elements
  return(gTree(children = gList(color_bar, border, ticks, labels, title), vp = legend_vp))
}

create_residual_legend_grob = function(max_res = 6) {
  # Define the same diverging color palette as the residual heatmap
  # Blue (#0571B0) for negative, White for zero, Red (#CA0020) for positive
  div_pal = colorRampPalette(c("#0571B0", "white", "#CA0020"))(100)
  
  # Create a matrix to draw the color gradient (rev to put Red at the top)
  grad_matrix = matrix(rev(div_pal), ncol = 1) 
  
  # Reduce height to 0.4 (40%) and shift left (x=0.2) to prevent cropping
  # Fix the width of the color bar itself to 5mm
  legend_vp = viewport(x = 0.2, y = 0.5, width = unit(5, "mm"), height = 0.4, just = "left")
  
  # Color bar body
  color_bar = rasterGrob(grad_matrix, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  border = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1))
  
  # Setup ticks and labels with 4 intervals (e.g., -6, -3, 0, 3, 6)
  y_ticks_at = seq(0, 1, by = 0.25) 
  labels_str = c(paste0("-", max_res), 
                 paste0("-", max_res / 2), 
                 "0", 
                 as.character(max_res / 2), 
                 as.character(max_res))
  
  ticks = segmentsGrob(
    x0 = unit(1, "npc"), y0 = unit(y_ticks_at, "npc"),
    x1 = unit(1, "npc") + unit(2, "mm"), y1 = unit(y_ticks_at, "npc"),
    gp = gpar(col = "black")
  )
  
  labels = textGrob(
    label = labels_str, 
    x = unit(1, "npc") + unit(3, "mm"), y = unit(y_ticks_at, "npc"), 
    just = "left", gp = gpar(fontsize = 10)
  )
  
  # Shift the title slightly upwards (+2.5 lines)
  title = textGrob("Residual", 
                   x = unit(0.5, "npc"), 
                   y = unit(1, "npc") + unit(2.5, "lines"), 
                   just = "center", gp = gpar(fontsize = 10, lineheight = 0.9))
  
  # Combine all elements
  return(gTree(children = gList(color_bar, border, ticks, labels, title), vp = legend_vp))
}
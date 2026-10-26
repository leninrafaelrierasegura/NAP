# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)

loglog_line_equation <- function(x1, y1, slope) {
  b <- log10(y1 / (x1 ^ slope))
  
  function(x) {
    (x ^ slope) * (10 ^ b)
  }
}
exp_line_equation <- function(x1, y1, slope) {
  lnC <- log(y1) - slope * x1
  
  function(x) {
    exp(lnC + slope * x)
  }
}
compute_guiding_lines <- function(x_axis_vector, errors, theoretical_rates, line_equation_fun) {
  guiding_lines <- matrix(NA, nrow = length(x_axis_vector), ncol = length(theoretical_rates))
  
  for (j in seq_along(theoretical_rates)) {
    guiding_lines_aux <- matrix(NA, nrow = length(x_axis_vector), ncol = length(x_axis_vector))
    
    for (k in seq_along(x_axis_vector)) {
      point_x1 <- x_axis_vector[k]
      point_y1 <- errors[k, j]
      slope <- theoretical_rates[j]
      
      line <- line_equation_fun(x1 = point_x1, y1 = point_y1, slope = slope)
      guiding_lines_aux[, k] <- line(x_axis_vector)
    }
    
    guiding_lines[, j] <- rowMeans(guiding_lines_aux)
  }
  
  return(guiding_lines)
}

# Function to order the vertices for plotting
plotting.order <- function(v, graph){
  edge_number <- graph$mesh$VtE[, 1]
  pos <- sum(edge_number == 1)+1
  return(c(v[1], v[3:pos], v[2], v[(pos+1):length(v)], v[2]))
}

# Function to set the scene for 3D plots
global.scene.setter <- function(x_range, y_range, z_range, z_aspectratio = 4) {
  
  return(list(xaxis = list(title = "x", range = x_range),
              yaxis = list(title = "y", range = y_range),
              zaxis = list(title = "z", range = z_range),
              aspectratio = list(x = 2*(1+2/pi), 
                                 y = 2*(2/pi), 
                                 z = z_aspectratio*(2/pi)),
              camera = list(eye = list(x = (1+2/pi)/2, 
                                       y = 4, 
                                       z = 2),
                            center = list(x = (1+2/pi)/2, 
                                          y = 0, 
                                          z = 0))))
}

global.scene.setter.rec.and.sphere <- function(x_range, y_range, z_range) {
  
  return(list(xaxis = list(title = list(text = "x", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = x_range),
              yaxis = list(title = list(text = "y", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = y_range),
              zaxis = list(title = list(text = "z", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = z_range),
              aspectratio = list(x = 1, 
                                 y = 1, 
                                 z = 1),
              camera = list(eye = list(x = 2, 
                                       y = 2, 
                                       z = 2),
                            center = list(x = 0, 
                                          y = 0, 
                                          z = 0))))
}


global.scene.setter.rec.and.sphere.for.hat <- function(x_range, y_range, z_range) {
  
  return(list(xaxis = list(title = list(text = "x", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = x_range),
              yaxis = list(title = list(text = "y", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = y_range),
              zaxis = list(title = list(text = "z", font = list(color = "gray", size = 16)),  tickfont = list(color = "gray", size = 13), range = z_range),
              aspectratio = list(x = 1, 
                                 y = 1, 
                                 z = 0.4),
              camera = list(eye = list(x = 2, 
                                       y = 2, 
                                       z = 2),
                            center = list(x = 0, 
                                          y = 0, 
                                          z = 0))))
}

# Function to plot in 3D
graph.plotter.3d.old <- function(graph, time_seq, frame_val_to_display, ...) {
  U_list <- list(...)
  U_names <- sapply(substitute(list(...))[-1], deparse)

  # Spatial coordinates
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  weights <- graph$mesh$weights

  # Apply plotting.order to each U
  U_list <- lapply(U_list, function(U) apply(U, 2, plotting.order, graph = graph))
  n_vars <- length(U_list)

  # Create plot_data frame with time and position replicated
  n_time <- ncol(U_list[[1]])
  base_data <- data.frame(
    x = rep(x, times = n_time),
    y = rep(y, times = n_time),
    the_graph = 0,
    frame = rep(time_seq, each = length(x))
  )

  # Add U columns to plot_data
  for (i in seq_along(U_list)) {
    base_data[[paste0("u", i)]] <- as.vector(U_list[[i]])
  }

  plot_data <- base_data

  # Generate vertical lines
  vertical_lines_list <- lapply(seq_along(U_list), function(i) {
    do.call(rbind, lapply(time_seq, function(t) {
      idx <- which(plot_data$frame == t)
      z_vals <- plot_data[[paste0("u", i)]][idx]
      data.frame(
        x = rep(plot_data$x[idx], each = 3),
        y = rep(plot_data$y[idx], each = 3),
        z = as.vector(t(cbind(0, z_vals, NA))),
        frame = rep(t, each = length(idx) * 3)
      )
    }))
  })

  # Set axis ranges
  z_range <- range(unlist(U_list))
  x_range <- range(x)
  y_range <- range(y)

  # Create plot
  p <- plot_ly(plot_data, frame = ~frame) %>%
    add_trace(x = ~x, y = ~y, z = ~the_graph, type = "scatter3d", mode = "lines",
              name = "", showlegend = FALSE,
              line = list(color = "black", width = 3))

  # Add traces for each variable
  colors <- rev(viridisLite::viridis(n_vars)) #RColorBrewer::brewer.pal(min(n_vars, 8), "Set1")
  for (i in seq_along(U_list)) {
    p <- add_trace(p,
      x = ~x, y = ~y, z = as.formula(paste0("~u", i)),
      type = "scatter3d", mode = "lines", name = U_names[i],
      line = list(color = colors[i], width = 3))
  }

  # Add vertical lines
  for (i in seq_along(vertical_lines_list)) {
    p <- add_trace(p,
      data = vertical_lines_list[[i]],
      x = ~x, y = ~y, z = ~z, frame = ~frame,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      name = "Vertical lines",
      showlegend = FALSE)
  }
  frame_name <- deparse(substitute(frame_val_to_display))
  # Layout and animation controls
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(type = "buttons", showactive = FALSE,
                              buttons = list(
                                list(label = "Play", method = "animate",
                                     args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE), fromcurrent = TRUE))),
                                list(label = "Pause", method = "animate",
                                     args = list(NULL, list(mode = "immediate", frame = list(duration = 0), redraw = FALSE)))
                              )
      )),
      title = paste0(frame_name,": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()

  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(title = paste0(frame_name,": ", formatC(frame_val_to_display[i], format = "f", digits = 4)))
  }

  return(p)
}

graph.plotter.3d <- function(graph, time_seq, frame_val_to_display, U_list) {
  U_names <- names(U_list) 
  # Spatial coordinates
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  weights <- graph$mesh$weights

  # Apply plotting.order to each U
  U_list <- lapply(U_list, function(U) apply(U, 2, plotting.order, graph = graph))
  n_vars <- length(U_list)
  
  # Create plot_data frame with time and position replicated
  n_time <- ncol(U_list[[1]])
  base_data <- data.frame(
    x = rep(x, times = n_time),
    y = rep(y, times = n_time),
    the_graph = 0,
    frame = rep(time_seq, each = length(x))
  )

  # Add U columns to plot_data
  for (i in seq_along(U_list)) {
    base_data[[paste0("u", i)]] <- as.vector(U_list[[i]])
  }

  plot_data <- base_data

  # Generate vertical lines
  vertical_lines_list <- lapply(seq_along(U_list), function(i) {
    do.call(rbind, lapply(time_seq, function(t) {
      idx <- which(plot_data$frame == t)
      z_vals <- plot_data[[paste0("u", i)]][idx]
      data.frame(
        x = rep(plot_data$x[idx], each = 3),
        y = rep(plot_data$y[idx], each = 3),
        z = as.vector(t(cbind(0, z_vals, NA))),
        frame = rep(t, each = length(idx) * 3)
      )
    }))
  })

  # Set axis ranges
  z_range <- range(unlist(U_list))
  x_range <- range(x)
  y_range <- range(y)

  # Create plot
  p <- plot_ly(plot_data, frame = ~frame) %>%
    add_trace(x = ~x, y = ~y, z = ~the_graph, type = "scatter3d", mode = "lines",
              name = "", showlegend = FALSE,
              line = list(color = "black", width = 3))

  if (n_vars == 2) {
    colors <- RColorBrewer::brewer.pal(min(n_vars, 8), "Set1") 
    } else {
    colors <- rev(viridisLite::viridis(n_vars)) 
  }
  # RColorBrewer::brewer.pal(min(n_vars, 8), "Set1")
  for (i in seq_along(U_list)) {
    p <- add_trace(p,
      x = ~x, y = ~y, z = as.formula(paste0("~u", i)),
      type = "scatter3d", mode = "lines", name = U_names[i],
      line = list(color = colors[i], width = 3))
  }

  # Add vertical lines
  for (i in seq_along(vertical_lines_list)) {
    p <- add_trace(p,
      data = vertical_lines_list[[i]],
      x = ~x, y = ~y, z = ~z, frame = ~frame,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      name = "Vertical lines",
      showlegend = FALSE)
  }
  frame_name <- deparse(substitute(frame_val_to_display))
  # Layout and animation controls
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(type = "buttons", showactive = FALSE,
                              buttons = list(
                                list(label = "Play", method = "animate",
                                     args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE), fromcurrent = TRUE))),
                                list(label = "Pause", method = "animate",
                                     args = list(NULL, list(mode = "immediate", frame = list(duration = 0), redraw = FALSE)))
                              )
      )),
      title = paste0(frame_name,": ", formatC(frame_val_to_display[1], format = "f", digits = 8))
    ) %>%
    plotly_build()

  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(title = paste0(frame_name,": ", formatC(frame_val_to_display[i], format = "f", digits = 8)))
  }

  return(p)
}

# Function to plot the error at each time step
error.at.each.time.plotter <- function(graph, U_true, U_approx, time_seq, time_step) {
  weights <- graph$mesh$weights
  error_at_each_time <- t(weights) %*% (U_true - U_approx)^2
  error <- sqrt(as.double(t(weights) %*% (U_true - U_approx)^2 %*% rep(time_step, ncol(U_true))))
  p <- plot_ly() %>% 
  add_trace(
  x = ~time_seq, y = ~error_at_each_time, type = 'scatter', mode = 'lines+markers',
  line = list(color = 'blue', width = 2),
  marker = list(color = 'blue', size = 4),
  name = "",
  showlegend = TRUE
) %>% 
  layout(
  title = paste0("Error at Each Time Step (Total error = ", formatC(error, format = "f", digits = 9), ")"),
  xaxis = list(title = "t"),
  yaxis = list(title = "Error"),
  legend = list(x = 0.1, y = 0.9)
)
  return(p)
}

# Function to plot the 3D comparison of U_true and U_approx
graph.plotter.3d.comparer <- function(graph, U_true, U_approx, time_seq) {
  x <- graph$mesh$V[, 1]; y <- graph$mesh$V[, 2]
  x <- plotting.order(x, graph); y <- plotting.order(y, graph)

  U_true <- apply(U_true, 2, plotting.order, graph = graph)
  U_approx <- apply(U_approx, 2, plotting.order, graph = graph)
  n_times <- length(time_seq)
  
  x_range <- range(x); y_range <- range(y); z_range <- range(c(U_true, U_approx))
  
  # Normalize time_seq
  time_normalized <- (time_seq - min(time_seq)) / (max(time_seq) - min(time_seq))
  blues <- colorRampPalette(c("lightblue", "blue"))(n_times)
  reds <- colorRampPalette(c("mistyrose", "red"))(n_times)
  
  # Accurate colorscales
  colorscale_greens <- Map(function(t, col) list(t, col), time_normalized, blues)
  colorscale_reds <- Map(function(t, col) list(t, col), time_normalized, reds)
  
  p <- plot_ly()
  
  # Static black graph structure
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 4),
              name = "Graph", showlegend = FALSE)
  
  # U_true traces (green)
  for (i in seq_len(n_times)) {
    z <- U_true[, i]
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x, y = y, z = z,
      line = list(color = blues[i], width = 4),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  
  # U_approx traces (dashed red)
  for (i in seq_len(n_times)) {
    z <- U_approx[, i]
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x, y = y, z = z,
      line = list(color = reds[i], width = 4, dash = "dot"),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  
  # Dummy green colorbar (True) – with ticks
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_greens,
    colorbar = list(
      title = list(font = list(size = 12, color = "black"), text = "Time", side = "top"),
      len = 0.9,
      thickness = 15,
      x = 1.02,
      xanchor = "left",
      y = 0.5,
      yanchor = "middle",
      tickvals = NULL,   # hide tick values
      ticktext = NULL,
      ticks = ""         # also hides tick marks
    ),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )

# Dummy red colorbar (Approx) – no ticks
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_reds,
    colorbar = list(
      title = list(font = list(size = 12, color = "black"), text = ".", side = "top"),
      len = 0.9,
      thickness = 15,
      x = 1.05,
      xanchor = "left",
      y = 0.5,
      yanchor = "middle"
    ),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 4),
              name = "Graph", showlegend = FALSE)
  p <- layout(p,
            scene = global.scene.setter(x_range, y_range, z_range),
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
  list(
    text = "Exact",
    x = 1.045,
    y = 0.5,
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    font = list(size = 12, color = "black"),
    textangle = -90
  ),
  list(
    text = "Approx",
    x = 1.075,
    y = 0.5,
    xref = "paper",
    yref = "paper",
    showarrow = FALSE,
    font = list(size = 12, color = "black"),
    textangle = -90
  )
)

)

  
  return(p)
}

# Function to plot a single 3D line for 
graph.plotter.3d.single <- function(graph, U_true, time_seq) {
  x <- graph$mesh$V[, 1]; y <- graph$mesh$V[, 2]
  x <- plotting.order(x, graph); y <- plotting.order(y, graph)

  U_true <- apply(U_true, 2, plotting.order, graph = graph)
  n_times <- length(time_seq)
  
  x_range <- range(x); y_range <- range(y); z_range <- range(U_true)
  z_range[1] <- z_range[1] - 10^-6
  viridis_colors <- viridisLite::viridis(100)
  
  # Normalize time_seq
  time_normalized <- (time_seq - min(time_seq)) / (max(time_seq) - min(time_seq))
  #greens <- colorRampPalette(c("palegreen", "darkgreen"))(n_times)
  greens <- colorRampPalette(c(viridis_colors[1], viridis_colors[50],  viridis_colors[100]))(n_times)
  # Accurate colorscales
  colorscale_greens <- Map(function(t, col) list(t, col), time_normalized, greens)
  
  p <- plot_ly()
  
  # Add the 3D lines with fading green color
  for (i in seq_len(n_times)) {
    z <- U_true[, i]
    
    p <- add_trace(
      p,
      type = "scatter3d",
      mode = "lines",
      x = x,
      y = y,
      z = z,
      line = list(color = greens[i], width = 2),
      showlegend = FALSE,
      scene = "scene"
    )
  }
  p <- p %>%
    add_trace(x = x, y = y, z = rep(0, length(x)),
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 5),
              name = "Graph", showlegend = FALSE)
  # Add dummy heatmap to show colorbar (not part of scene)
  p <- add_trace(
    p,
    type = "heatmap",
    z = matrix(time_seq, nrow = 1),
    showscale = TRUE,
    colorscale = colorscale_greens,
    colorbar = list(
    title = list(font = list(size = 12, color = "black"), text = "Time", side = "top"),
    len = 0.9,         # height (0 to 1)
    thickness = 15,     # width in pixels
    x = 1.02,           # shift it slightly right of the plot
    xanchor = "left",
    y = 0.5,
    yanchor = "middle"),
    x = matrix(time_seq, nrow = 1),
    y = matrix(1, nrow = 1),
    hoverinfo = "skip",
    opacity = 0
  )
  
  p <- layout(p,
              scene = global.scene.setter(x_range, y_range, z_range),
              xaxis = list(visible = FALSE),
              yaxis = list(visible = FALSE)
  )
  
  return(p)
}

# Function to plot the error convergence
error.convergence.plotter <- function(x_axis_vector, 
                                      alpha_vector, 
                                      errors, 
                                      theoretical_rates, 
                                      observed_rates,
                                      line_equation_fun,
                                      fig_title,
                                      x_axis_label,
                                      apply_sqrt = FALSE) {
  
  x_vec <- if (apply_sqrt) sqrt(x_axis_vector) else x_axis_vector
  
  guiding_lines <- compute_guiding_lines(x_axis_vector = x_vec, 
                                         errors = errors, 
                                         theoretical_rates = theoretical_rates, 
                                         line_equation_fun = line_equation_fun)
  
  default_colors <- scales::hue_pal()(length(alpha_vector))
  
  plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
    geom_line(
      data = data.frame(x = x_vec, y = guiding_lines[, i]),
      aes(x = x, y = y),
      color = default_colors[i],
      linetype = "dashed",
      show.legend = FALSE
    )
  })
  
  df <- as.data.frame(cbind(x_vec, errors))
  colnames(df) <- c("x_axis_vector", alpha_vector)
  df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
  
  custom_labels <- paste0(formatC(alpha_vector, format = "f", digits = 2), 
                          " | ", 
                          formatC(theoretical_rates, format = "f", digits = 4), 
                          " | ", 
                          formatC(observed_rates, format = "f", digits = 4))
  
  df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)

  p <- ggplot() +
    geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    plot_lines +
    labs(
      title = fig_title,
      x = x_axis_label,
      y = expression(Error),
      color = if (apply_sqrt) "alpha | bound benchmark | observed" else "alpha | theoretical | observed"
    ) +
    (if (apply_sqrt) {
      scale_x_continuous(breaks = x_vec, labels = round(x_axis_vector, 4))
    } else {
      scale_x_log10(breaks = x_axis_vector, labels = round(x_axis_vector, 4))
    }) +
    (if (apply_sqrt) {
      scale_y_continuous(trans = "log", labels = scales::scientific_format())
    } else {
      scale_y_log10(labels = scales::scientific_format())
    }) +
    theme_minimal() +
    theme(text = element_text(family = "Palatino"),
          legend.position = "bottom",
          legend.direction = "vertical",
          plot.margin = margin(0, 0, 0, 0),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  
  return(p)
}

errorConvergencePlotterLatex <- function(x_axis_vector, 
                                      alpha_vector, 
                                      errors, 
                                      theoretical_rates, 
                                      observed_rates,
                                      line_equation_fun,
                                      fig_title,
                                      x_axis_label,
                                      y_axis_label,
                                      color_label,
                                      apply_sqrt = FALSE) {

  relative_per_error <- 100 * abs(theoretical_rates - observed_rates) / abs(theoretical_rates)
  
  x_vec <- if (apply_sqrt) sqrt(x_axis_vector) else x_axis_vector
  
  guiding_lines <- compute_guiding_lines(x_axis_vector = x_vec, 
                                         errors = errors, 
                                         theoretical_rates = theoretical_rates, 
                                         line_equation_fun = line_equation_fun)
  
  default_colors <- scales::hue_pal()(length(alpha_vector))
  
  plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
    geom_line(
      data = data.frame(x = x_vec, y = guiding_lines[, i]),
      aes(x = x, y = y),
      color = default_colors[i],
      linetype = "dashed",
      linewidth = 0.6,
      show.legend = FALSE
    )
  })
  
  df <- as.data.frame(cbind(x_vec, errors))
  colnames(df) <- c("x_axis_vector", alpha_vector)
  df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
  if (apply_sqrt) {
    custom_labels <- sprintf("$%s \\; | \\; %s \\; | \\; %s$",
      formatC(alpha_vector, format = "f", digits = 2),
      formatC(theoretical_rates, format = "f", digits = 4),
      formatC(observed_rates, format = "f", digits = 4))
  } else {
    custom_labels <- sprintf("$%s \\; | \\; %s \\; | \\; %s \\; | \\; %s$",
      formatC(alpha_vector, format = "f", digits = 2),
      formatC(theoretical_rates, format = "f", digits = 4),
      formatC(observed_rates, format = "f", digits = 4),
      formatC(relative_per_error, format = "f", digits = 2))
  }

  df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)

  scientific_ticks <- function(x) {
    labels <- formatC(x, format = "e", digits = 2)
    labels <- sub("0+e", "e", labels)
    labels <- sub("[.]e", "e", labels)
    labels <- sub("e([+-])0+([0-9]+)$", "e\\1\\2", labels)
    sub("e[+]", "e", labels)
  }

  p <- ggplot() +
    geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column), linewidth = 1.6) +
    geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column), size = 2) +
    plot_lines +
    labs(
      title = fig_title,
      x = x_axis_label,
      y = y_axis_label,
      color = color_label
    ) +
    (if (apply_sqrt) {
      scale_x_continuous(breaks = x_vec, labels = formatC(x_axis_vector, format = "f", digits = 0))
    } else {
      scale_x_log10(breaks = x_axis_vector, labels = scientific_ticks)
    }) +
    (if (apply_sqrt) {
      #scale_y_continuous(trans = "log", labels = scales::scientific_format())
      scale_y_continuous(trans = scales::log_trans(base = exp(1)), labels = scales::scientific_format())
    } else {
      scale_y_log10(labels = scales::scientific_format())
    }) +
    theme_minimal() +
    # theme(text = element_text(family = "Palatino"),
    #       legend.position = "bottom",
    #       legend.direction = "vertical",
    #       #plot.margin = margin(0, 0, 0, 0),
    #       plot.title = element_text(hjust = 0.5, size = 18, face = "plain")) +
    theme(
  text = element_text(family = "Palatino", size = 18),

  plot.title = element_text(
    hjust = 0.5,
    size = 20,
    face = "plain"
  ),

  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),

  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),

  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),

  legend.position = "bottom",
  legend.direction = "vertical"
)
  
  return(p)
}

graph.plotter.3d.static <- function(graph, z_list) {
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  U_names <- names(z_list)
  n_vars <- length(z_list)
  z_list <- lapply(z_list, function(z) plotting.order(z, graph))

  # Axis ranges
  z_range <- range(unlist(z_list))
  x_range <- range(x)
  y_range <- range(y)

  if (n_vars == 2) {
    colors <- RColorBrewer::brewer.pal(min(n_vars, 8), "Set1") 
    } else {
    colors <- rev(viridisLite::viridis(n_vars)) 
  }
  p <- plot_ly()

  for (i in seq_along(z_list)) {
    z <- z_list[[i]]

    # Main 3D curve
    p <- add_trace(
      p,
      x = x, y = y, z = z,
      type = "scatter3d", mode = "lines",
      line = list(color = colors[i], width = 3),
      name = U_names[i], showlegend = TRUE
    )

    # Efficient vertical lines: one trace with breaks (NA)
    x_vert <- rep(x, each = 3)
    y_vert <- rep(y, each = 3)
    z_vert <- unlist(lapply(z, function(zj) c(0, zj, NA)))

    p <- add_trace(
      p,
      x = x_vert, y = y_vert, z = z_vert,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      showlegend = FALSE
    )
  }
  p <- p %>% add_trace(x = x, y = y, z = x*0, type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 3),
              name = "thegraph", showlegend = FALSE) %>%
    layout(scene = global.scene.setter(x_range, y_range, z_range))
  return(p)
}


graph.plotter.3d.two.meshes.time <- function(graph_finer, graph_coarser, 
                                             time_seq, frame_val_to_display,
                                             fs_finer = list(), fs_coarser = list()) {
  # Spatial coordinates (ordered for plotting)
  x_finer <- plotting.order(graph_finer$mesh$V[, 1], graph_finer)
  y_finer <- plotting.order(graph_finer$mesh$V[, 2], graph_finer)
  x_coarser <- plotting.order(graph_coarser$mesh$V[, 1], graph_coarser)
  y_coarser <- plotting.order(graph_coarser$mesh$V[, 2], graph_coarser)
  
  n_time <- if (length(fs_finer) > 0) ncol(fs_finer[[1]]) else ncol(fs_coarser[[1]])

  # Helper: make dataframe from one function
  make_df <- function(f_mat, graph, x, y, mesh_name) {
    z <- apply(f_mat, 2, plotting.order, graph = graph)
    data.frame(
      x = rep(x, times = n_time),
      y = rep(y, times = n_time),
      z = as.vector(z),
      frame = rep(time_seq, each = length(x)),
      mesh = mesh_name
    )
  }
  
  # Build data for finer functions
  data_finer_list <- lapply(names(fs_finer), function(nm) {
    make_df(fs_finer[[nm]], graph_finer, x_finer, y_finer, nm)
  })
  
  # Build data for coarser functions
  data_coarser_list <- lapply(names(fs_coarser), function(nm) {
    make_df(fs_coarser[[nm]], graph_coarser, x_coarser, y_coarser, nm)
  })
  
  # Combine
  all_data <- c(data_finer_list, data_coarser_list)
  
  # Baseline graph (on finer mesh for consistency)
  data_graph <- data.frame(
    x = rep(x_finer, times = n_time),
    y = rep(y_finer, times = n_time),
    z = 0,
    frame = rep(time_seq, each = length(x_finer)),
    mesh = "Graph"
  )
  
# --------- Vertical lines helper ----------
vertical_lines <- function(x, y, z, frame_vals, mesh_name) {
  do.call(rbind, lapply(seq_along(frame_vals), function(i) {
    idx <- ((i - 1) * length(x) + 1):(i * length(x))
    data.frame(
      x = rep(x, each = 3),
      y = rep(y, each = 3),
      z = as.vector(t(cbind(0, z[idx], NA))),
      frame = rep(frame_vals[i], each = length(x) * 3),
      mesh = mesh_name
    )
  }))
}

# --------- Compute vertical lines per mesh using max absolute value ---------
make_vertical_from_list <- function(data_list, x, y, mesh_name) {
  if (length(data_list) == 0) return(NULL)
  
  # Reshape each function's z back to matrix: (nodes × time)
  z_mats <- lapply(data_list, function(df) {
    matrix(df$z, nrow = length(x), ncol = length(time_seq))
  })
  
  # Stack into 3D array: (nodes × time × functions)
  arr <- array(unlist(z_mats), dim = c(length(x), length(time_seq), length(z_mats)))
  
  # For each node × time, select the entry with largest absolute value (keep sign)
  idx <- apply(arr, c(1, 2), function(v) which.max(abs(v)))
  z_signed_max <- mapply(function(i, j) arr[i, j, idx[i, j]],
                         rep(1:length(x), times = length(time_seq)),
                         rep(1:length(time_seq), each = length(x)))
  
  # Flatten back into long vector
  z_signed_max <- as.vector(z_signed_max)
  
  vertical_lines(x, y, z_signed_max, time_seq, mesh_name)
}


vertical_finer   <- make_vertical_from_list(data_finer_list,   x_finer,   y_finer,   "finer")
vertical_coarser <- make_vertical_from_list(data_coarser_list, x_coarser, y_coarser, "coarser")

  
  # Compute ranges
  all_z <- unlist(lapply(all_data, function(df) df$z))
  x_range <- range(c(x_finer, x_coarser))
  y_range <- range(c(y_finer, y_coarser))
  z_range <- range(all_z)
  
  # --------- Plotly object ----------
  p <- plot_ly(frame = ~frame)
  
  # Add traces for finer + coarser (looping automatically with names)
  for (df in all_data) {
    p <- p %>%
      add_trace(data = df,
                x = ~x, y = ~y, z = ~z,
                type = "scatter3d", mode = "lines",
                line = list(width = 3),
                name = unique(df$mesh))
  }
  
  # Add baseline
  p <- p %>%
    add_trace(data = data_graph,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "black", width = 2),
              name = "Graph", showlegend = FALSE)
  
# Add verticals (one per mesh, envelope of all functions)
if (!is.null(vertical_finer)) {
  p <- p %>%
    add_trace(data = vertical_finer,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "gray", width = 0.5),
              name = "Vertical finer", showlegend = FALSE)
}
if (!is.null(vertical_coarser)) {
  p <- p %>%
    add_trace(data = vertical_coarser,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d", mode = "lines",
              line = list(color = "gray", width = 0.5),
              name = "Vertical coarser", showlegend = FALSE)
}

  
  frame_name <- deparse(substitute(frame_val_to_display))
  
  p <- p %>%
    layout(
      scene = global.scene.setter(x_range, y_range, z_range),
      updatemenus = list(list(
        type = "buttons", showactive = FALSE,
        buttons = list(
          list(label = "Play", method = "animate",
               args = list(NULL, list(frame = list(duration = 2000 / length(time_seq), redraw = TRUE),
                                      fromcurrent = TRUE))),
          list(label = "Pause", method = "animate",
               args = list(NULL, list(mode = "immediate", 
                                      frame = list(duration = 0), redraw = FALSE)))
        )
      )),
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[1], format = "f", digits = 4))
    ) %>%
    plotly_build()
  
  # Update frame titles
  for (i in seq_along(p$x$frames)) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", formatC(frame_val_to_display[i], format = "f", digits = 4))
    )
  }
  
  return(p)
}

plot_3d_rectangle_scatter <- function(loc, eigvals, eigfuncs, fixed_colorscale = TRUE) {

  x <- loc[,1]
  y <- loc[,2]

  colorscale <- "Viridis"

  # Compute global color limits if fixed
  if (fixed_colorscale) {
    cmin <- min(eigfuncs)
    cmax <- max(eigfuncs)
  } else {
    cmin <- NULL
    cmax <- NULL
  }

  # Frames ----------------------------------------------------------------
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = eigfuncs[,i],
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 5,
          color = eigfuncs[,i],
          colorscale = colorscale,
          showscale = TRUE,
          cmin = cmin,
          cmax = cmax
        )
      ))
    )
  })

  # Initial Plot ----------------------------------------------------------
  p <- plot_ly(
    x = x,
    y = y,
    z = eigfuncs[,1],
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = eigfuncs[,1],
      colorscale = colorscale,
      showscale = TRUE,
      cmin = cmin,
      cmax = cmax
    ),
    frame = "1"
  )

  p$x$frames <- frames

  z_range <- range(eigfuncs)
  x_range <- range(x)
  y_range <- range(y)

  frame_name <- deparse(substitute(eigvals))

  # Layout + slider -------------------------------------------------------
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(
              list(as.character(i)),
              list(frame = list(duration = 300, redraw = TRUE),
                   mode = "immediate")
            )
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(
                 NULL,
                 list(frame = list(duration = 300, redraw = TRUE),
                      fromcurrent = TRUE)
               )),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0))))
        )
      )
    )
  ) %>% plotly_build()

  # Update title in each frame --------------------------------------------
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }

  return(p)
}

plot_3d_rectangle_onecol_vir <- function(mesh, fvals) {
  
  colorscale = "Viridis"
  opacity = 1
  
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  tri <- mesh$graph$tv
  
  # fvals is a vector
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(fvals)
  
  # Static mesh3d plot
  p <- plot_ly(
    x = x, y = y, z = fvals,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    intensity = fvals,
    colorscale = colorscale,
    opacity = opacity,
    flatshading = TRUE,
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "f: ", sprintf("%.5f", fvals)
    ),
    hoverinfo = "text"
  ) %>% layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}

plot_3d_rectangle_onecol <- function(mesh, fvals) {
  
  opacity <- 1
  surface_color <- "blue"
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  tri <- mesh$graph$tv
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(fvals)
  
  # Static mesh3d plot
  p <- plot_ly(
    x = x, y = y, z = fvals,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    color = surface_color,     # Fixed blue color
    opacity = opacity,
    flatshading = TRUE,
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "f: ", sprintf("%.5f", fvals)
    ),
    hoverinfo = "text"
  ) %>% layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}

plot_3d_square_mesh <- function(mesh) {
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  # Flat square (or use mesh$loc[,3] if available)
  if (ncol(mesh$loc) >= 3) {
    z <- mesh$loc[,3]
  } else {
    z <- rep(0, length(x))
  }
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- c(0,1)
  
  tri <- mesh$graph$tv
  
  # ---- Uniform gray color for faces ----
  gray_rgb <- rep(list(rgb(180, 180, 180, maxColorValue = 255)), length(x))
  white_rgb <- rep(list(rgb(255, 255, 255, alpha = 255, maxColorValue = 255)), length(x))

  
  df3 <- data.frame(x = x, 
                  y = y, 
                  z = z)
  
  # ---- Plot mesh faces ----
  p <- plot_ly() %>% add_trace(
    x = x,
    y = y,
    z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    vertexcolor = white_rgb,
    flatshading = TRUE,
    opacity = 1,
    showscale = FALSE
  )
  
  # ---- Extract unique edges ----
  edges <- rbind(
    tri[,1:2],
    tri[,2:3],
    tri[,c(3,1)]
  )
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  # ---- Efficient edges as one trace using NA breaks ----
  x_edges <- as.vector(t(cbind(x[edges[,1]], x[edges[,2]], NA)))
  y_edges <- as.vector(t(cbind(y[edges[,1]], y[edges[,2]], NA)))
  z_edges <- as.vector(t(cbind(z[edges[,1]], z[edges[,2]], NA)))
  
  p <- add_trace(
    p,
    x = x_edges, y = y_edges, z = z_edges,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 3),
    showlegend = FALSE,
    hoverinfo = "none"
  ) %>% 
    add_trace(data = df3, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = "gray", symbol = 104)) %>%
    layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}


plot_3d_square_mesh_with_hat <- function(mesh,
                                hat_nodes = NULL,
                                hat_height = 1,
                                hat_color = NULL,
                                hat_alpha = 0.6) {
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  n <- length(x)
  
  z0 <- rep(0, n)
  tri <- mesh$graph$tv
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- c(0, hat_height)
  
  # ---- Extract unique edges once ----
  edges <- rbind(tri[,1:2], tri[,2:3], tri[,c(3,1)])
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  # Base mesh color
  white_rgb <- rep(list(rgb(255, 255, 255, alpha = 255, maxColorValue = 255)), n)
  
  df3 <- data.frame(x = x, 
                  y = y, 
                  z = rep(0, length(x)))
  
  # ---- Base mesh ----
  p <- plot_ly() %>% add_trace(
    x = x, y = y, z = z0,
    i = tri[,1] - 1, j = tri[,2] - 1, k = tri[,3] - 1,
    type = "mesh3d",
    vertexcolor = white_rgb,
    flatshading = TRUE,
    opacity = 1,
    name = "Mesh"
  )
  
  # ---- Base mesh edges ----
  x_edges <- as.vector(t(cbind(x[edges[,1]], x[edges[,2]], NA)))
  y_edges <- as.vector(t(cbind(y[edges[,1]], y[edges[,2]], NA)))
  z_edges <- as.vector(t(cbind(z0[edges[,1]], z0[edges[,2]], NA)))
  
  # ---- Hat functions ----
  if (!is.null(hat_nodes)) {
    
    if (is.null(hat_color)) {
      hat_color <- grDevices::rainbow(length(hat_nodes))
    }
    if (length(hat_color) != length(hat_nodes)) {
      stop("length(hat_color) must equal length(hat_nodes)")
    }
    
    for (k in seq_along(hat_nodes)) {
      i0 <- hat_nodes[k]
      
      phi <- rep(0, n)
      phi[i0] <- 1
      z_hat <- hat_height * phi
      
      hat_rgb <- rep(list(adjustcolor(hat_color[k], alpha.f = hat_alpha)), n)
      
      # ---- Hat surface ----
      p <- add_trace(p,
        x = x, y = y, z = z_hat,
        i = tri[,1] - 1, j = tri[,2] - 1, k = tri[,3] - 1,
        type = "mesh3d",
        vertexcolor = hat_rgb,
        flatshading = TRUE,
        opacity = 1,
        name = paste0("phi_", i0),
        showlegend = FALSE
      )
      
      # ---- Hat edges (wireframe) ----
      x_hat_edges <- as.vector(t(cbind(x[edges[,1]], x[edges[,2]], NA)))
      y_hat_edges <- as.vector(t(cbind(y[edges[,1]], y[edges[,2]], NA)))
      z_hat_edges <- as.vector(t(cbind(z_hat[edges[,1]], z_hat[edges[,2]], NA)))
      
      p <- add_trace(p,
        x = x_hat_edges, y = y_hat_edges, z = z_hat_edges,
        type = "scatter3d", mode = "lines",
        line = list(color = hat_color[k], width = 2),
        showlegend = FALSE
      )
      
      # # Mark hat center node
      # p <- add_trace(p,
      #   x = x[i0], y = y[i0], z = z_hat[i0],
      #   type = "scatter3d", mode = "markers",
      #   marker = list(size = 8, color = hat_color[k]),
      #   name = paste0("node_", i0)
      # )
    }
  }
  
  p <- add_trace(p,
    x = x_edges, y = y_edges, z = z_edges,
    type = "scatter3d", mode = "lines",
    line = list(color = "black", width = 4),
    showlegend = FALSE
  )
  
  p <- p %>% 
    add_trace(data = df3, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = "gray", symbol = 104), showlegend = FALSE) %>% 
    layout(font = list(family = "Palatino Linotype"),
    scene = global.scene.setter.rec.and.sphere.for.hat(x_range, y_range, z_range)
  )
  
  return(p)
}

plot_3d_sphere_scatter <- function(mesh, 
                                          eigvals, 
                                          eigfuncs,
                                          fixed_colorscale = TRUE) {
  
  colorscale = "Viridis"
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  
  # Compute global color limits if fixed
  if (fixed_colorscale) {
    cmin <- min(eigfuncs)
    cmax <- max(eigfuncs)
  } else {
    cmin <- NULL
    cmax <- NULL
  }
  
  # Create frames
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    
    fvals <- eigfuncs[,i]
    
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 5,
          color = fvals,
          colorscale = colorscale,
          showscale = TRUE,
          cmin = cmin,
          cmax = cmax
        ),
        text = paste0(
          "x: ", sprintf("%.3f", x), "<br>",
          "y: ", sprintf("%.3f", y), "<br>",
          "z: ", sprintf("%.3f", z), "<br>",
          "f: ", sprintf("%.5f", fvals)
        ),
        hoverinfo = "text"
      ))
    )
  })
  
  # Initial plot (frame 1)
  fvals0 <- eigfuncs[,1]
  
  p <- plot_ly(
    x = x, y = y, z = z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = fvals0,
      colorscale = colorscale,
      showscale = TRUE,
      cmin = cmin,
      cmax = cmax
    ),
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals0)
    ),
    hoverinfo = "text",
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(z)
  
  # Slider + play/pause
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Mode: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(
            label = "Play",
            method = "animate",
            args = list(NULL,
                        list(frame = list(duration = 300, redraw = TRUE),
                             fromcurrent = TRUE,
                             mode = "immediate"))
          ),
          list(
            label = "Pause",
            method = "animate",
            args = list(NULL,
                        list(frame = list(duration = 0, redraw = FALSE),
                             mode = "immediate"))
          )
        )
      )
    )
  ) %>% plotly_build()
  
  # Update title per frame
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}

plot_3d_slider_sphere <- function(mesh, eigvals, eigfuncs, fixed_colorscale = TRUE) {
  
  colorscale = "Viridis"
  opacity = 1
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  z <- mesh$loc[, 3]
  tri <- mesh$graph$tv
  
  # Global intensity limits if fixed
  if (fixed_colorscale) {
    cmin <- min(eigfuncs)
    cmax <- max(eigfuncs)
  } else {
    cmin <- NULL
    cmax <- NULL
  }
  
  # Create frames
  frames <- lapply(seq_len(ncol(eigfuncs)), function(i) {
    fvals <- eigfuncs[,i]
    list(
      name = as.character(i),
      data = list(list(
        x = x,
        y = y,
        z = z,
        i = tri[,1] - 1,
        j = tri[,2] - 1,
        k = tri[,3] - 1,
        type = "mesh3d",
        intensity = fvals,
        colorscale = colorscale,
        opacity = opacity,
        flatshading = TRUE,
        cmin = cmin,
        cmax = cmax,
        text = paste0(
          "x: ", sprintf("%.3f", x), "<br>",
          "y: ", sprintf("%.3f", y), "<br>",
          "z: ", sprintf("%.3f", z), "<br>",
          "f: ", sprintf("%.5f", fvals)
        ),
        hoverinfo = "text"
      ))
    )
  })
  
  # Initial plot
  fvals0 <- eigfuncs[,1]
  
  p <- plot_ly(
    x = x, y = y, z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    intensity = fvals0,
    colorscale = colorscale,
    opacity = opacity,
    flatshading = TRUE,
    cmin = cmin,
    cmax = cmax,
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals0)
    ),
    colorbar = list(
      thickness = 15,
      len = 0.5,      # fraction of the height of the plot
      y = 0.5,        # center vertically (0 = bottom, 1 = top)
      yanchor = "middle",
      title = ""
    ),
    hoverinfo = "text",
    frame = "1"
  )
  
  frame_name <- deparse(substitute(eigvals))
  
  p$x$frames <- frames
  
  # Layout + slider + buttons
  p <- p %>% layout(
    title = paste0(frame_name, ": ", eigvals[1]),
    scene = global.scene.setter.rec.and.sphere(range(x), range(y), range(z)),
    sliders = list(
      list(
        active = 0,
        currentvalue = list(prefix = "Frame: "),
        pad = list(t = 50),
        steps = lapply(seq_len(ncol(eigfuncs)), function(i) {
          list(
            label = as.character(i),
            method = "animate",
            args = list(list(as.character(i)),
                        list(mode = "immediate",
                             frame = list(duration = 300, redraw = TRUE),
                             transition = list(duration = 0)))
          )
        })
      )
    ),
    updatemenus = list(
      list(
        type = "buttons",
        showactive = FALSE,
        y = 1,
        x = 1.15,
        xanchor = "right",
        yanchor = "top",
        buttons = list(
          list(label = "Play",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 300, redraw = TRUE),
                                      fromcurrent = TRUE, mode = "immediate"))),
          list(label = "Pause",
               method = "animate",
               args = list(NULL, list(frame = list(duration = 0, redraw = FALSE),
                                      mode = "immediate")))
        )
      )
    )
  ) %>% plotly_build()
  
  # Add titles for frames
  for (i in seq_len(ncol(eigfuncs))) {
    p$x$frames[[i]]$layout <- list(
      title = paste0(frame_name, ": ", eigvals[i])
    )
  }
  
  return(p)
}


plot_3d_sphere_scatter_onecol <- function(mesh, eigfuncs) {
  
  colorscale = "Viridis"
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  
  # eigfuncs is a vector
  fvals <- eigfuncs
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(z)
  
  # Build static 3D scatter plot
  p <- plot_ly(
    x = x, y = y, z = z,
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 5,
      color = fvals,
      colorscale = colorscale,
      showscale = TRUE
    ),
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals)
    ),
    hoverinfo = "text"
  ) %>% layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}



# plot_3d_sphere_surface_onecol <- function(mesh, fvals) {
#   
#   colorscale = "Viridis"
#   opacity = 1
#   
#   
#   x <- mesh$loc[, 1]
#   y <- mesh$loc[, 2]
#   z <- mesh$loc[, 3]
#   tri <- mesh$graph$tv
#   
#   # fvals is a vector
#   
#   x_range <- range(x)
#   y_range <- range(y)
#   z_range <- range(z)
#   
#   # Static mesh3d plot
#   p <- plot_ly(
#     x = x, y = y, z = z,
#     i = tri[,1] - 1,
#     j = tri[,2] - 1,
#     k = tri[,3] - 1,
#     type = "mesh3d",
#     intensity = fvals,
#     colorscale = colorscale,
#     opacity = opacity,
#     flatshading = TRUE,
#     text = paste0(
#       "x: ", sprintf("%.3f", x), "<br>",
#       "y: ", sprintf("%.3f", y), "<br>",
#       "z: ", sprintf("%.3f", z), "<br>",
#       "f: ", sprintf("%.5f", fvals)
#     ),
#     hoverinfo = "text"
#   ) %>% layout(
#     scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
#   )
#   
#   return(p)
# }


plot_3d_sphere_surface_onecol <- function(mesh, fvals) {
  
  colorscale = "Viridis"
  opacity = 1
  
  x <- mesh$loc[, 1]
  y <- mesh$loc[, 2]
  z <- mesh$loc[, 3]
  tri <- mesh$graph$tv
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(z)
  
  # Static mesh3d plot
  p <- plot_ly(
    x = x, y = y, z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    intensity = fvals,
    colorscale = colorscale,
    opacity = opacity,
    flatshading = TRUE,
    text = paste0(
      "x: ", sprintf("%.3f", x), "<br>",
      "y: ", sprintf("%.3f", y), "<br>",
      "z: ", sprintf("%.3f", z), "<br>",
      "f: ", sprintf("%.5f", fvals)
    ),
    hoverinfo = "text",
    colorbar = list(
      thickness = 15,
      len = 0.5,      # fraction of the height of the plot
      y = 0.5,        # center vertically (0 = bottom, 1 = top)
      yanchor = "middle",
      title = ""
    )
  ) %>% layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}


plot_3d_mesh_edges_faces_old <- function(mesh) {
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  tri <- mesh$graph$tv
  
  # ---- Uniform gray color for faces (vertexcolor) ----
  gray_rgb <- rep(list(rgb(180, 180, 180, maxColorValue = 255)), length(x))
  
  p <- plot_ly() %>% add_trace(
    x = x,
    y = y,
    z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    vertexcolor = gray_rgb,   # <- the correct way
    flatshading = TRUE,
    opacity = 1,
    showscale = FALSE
  )
  
  # ---- Extract unique edges ----
  edges <- rbind(
    tri[,1:2],
    tri[,2:3],
    tri[,c(3,1)]
  )
  
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  # ---- Add edges as blue lines ----
  for (e in 1:nrow(edges)) {
    i1 <- edges[e,1]
    i2 <- edges[e,2]
    
    p <- p %>% add_trace(
      x = c(x[i1], x[i2]),
      y = c(y[i1], y[i2]),
      z = c(z[i1], z[i2]),
      type = "scatter3d",
      mode = "lines",
      line = list(color = "blue", width = 3),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  
  return(p)
}

plot_3d_mesh_edges_faces <- function(mesh) {
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  tri <- mesh$graph$tv
  
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(z)
  
  # ---- Uniform gray color for faces ----
  gray_rgb <- rep(list(rgb(180, 180, 180, maxColorValue = 255)), length(x))
  white_rgb <- rep(list(rgb(255, 255, 255, alpha = 255, maxColorValue = 255)), length(x))
  
  df3 <- data.frame(x = x, 
                  y = y, 
                  z = z)
  
  # ---- Plot mesh faces ----
  p <- plot_ly() %>% add_trace(
    x = x,
    y = y,
    z = z,
    i = tri[,1] - 1,
    j = tri[,2] - 1,
    k = tri[,3] - 1,
    type = "mesh3d",
    vertexcolor = white_rgb,
    flatshading = TRUE,
    opacity = 1,
    showscale = FALSE
  )
  
  # ---- Extract unique edges ----
  edges <- rbind(
    tri[,1:2],
    tri[,2:3],
    tri[,c(3,1)]
  )
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  # ---- Efficient edges as a single trace using NA breaks ----
  x_edges <- as.vector(t(cbind(x[edges[,1]], x[edges[,2]], NA)))
  y_edges <- as.vector(t(cbind(y[edges[,1]], y[edges[,2]], NA)))
  z_edges <- as.vector(t(cbind(z[edges[,1]], z[edges[,2]], NA)))
  
  p <- add_trace(
    p,
    x = x_edges, y = y_edges, z = z_edges,
    type = "scatter3d",
    mode = "lines",
    line = list(color = "black", width = 3),
    showlegend = FALSE,
    hoverinfo = "none"
  ) %>% 
    add_trace(data = df3, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = "gray", symbol = 104)) %>%
    layout(
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}

plot_3d_mesh_edges_faces_with_hat <- function(mesh,
                                              hat_nodes = NULL,
                                              hat_height = 0.1,
                                              hat_color = NULL,
                                              hat_alpha = 0.6) {
  
  x <- mesh$loc[,1]
  y <- mesh$loc[,2]
  z <- mesh$loc[,3]
  n <- length(x)
  tri <- mesh$graph$tv
  
  x_range <- range(x)*2
  y_range <- range(y)*2
  z_range <- range(z)*2
  
  # ---- Extract unique edges once ----
  edges <- rbind(tri[,1:2], tri[,2:3], tri[,c(3,1)])
  edges <- t(apply(edges, 1, sort))
  edges <- unique(edges)
  
  # Base mesh colors
  white_rgb <- rep(list(rgb(255, 255, 255, alpha = 255, maxColorValue = 255)), n)
  
  df3 <- data.frame(x = x, y = y, z = z)
  
  # ---- Base spherical mesh ----
  p <- plot_ly() %>% add_trace(
    x = x, y = y, z = z,
    i = tri[,1] - 1, j = tri[,2] - 1, k = tri[,3] - 1,
    type = "mesh3d",
    vertexcolor = white_rgb,
    flatshading = TRUE,
    opacity = 1,
    showscale = FALSE,
    name = "Sphere mesh"
  )
  
  # ---- Base mesh edges (single efficient trace) ----
  x_edges <- as.vector(t(cbind(x[edges[,1]], x[edges[,2]], NA)))
  y_edges <- as.vector(t(cbind(y[edges[,1]], y[edges[,2]], NA)))
  z_edges <- as.vector(t(cbind(z[edges[,1]], z[edges[,2]], NA)))
  
 
  
  # ---- Hat functions on sphere ----
  if (!is.null(hat_nodes)) {
    
    # Validate indices
    if (any(hat_nodes < 1 | hat_nodes > n)) {
      stop("hat_nodes contains invalid node indices")
    }
    
    # Color recycling
    if (is.null(hat_color)) {
      hat_color <- grDevices::rainbow(length(hat_nodes))
    }
    hat_color <- rep(hat_color, length.out = length(hat_nodes))
    
    for (k in seq_along(hat_nodes)) {
      i0 <- hat_nodes[k]
      
      # Nodal hat basis
      phi <- rep(0, n)
      phi[i0] <- 1
      
      # Lift radially outward (important for sphere!)
      r <- sqrt(x^2 + y^2 + z^2)
      nx <- x / r
      ny <- y / r
      nz <- z / r
      
      x_hat <- x + hat_height * phi * nx
      y_hat <- y + hat_height * phi * ny
      z_hat <- z + hat_height * phi * nz
      
      # Hat surface color
      hat_rgb <- rep(list(adjustcolor(hat_color[k], alpha.f = hat_alpha)), n)
      
      # ---- Hat surface ----
      p <- add_trace(p,
        x = x_hat, y = y_hat, z = z_hat,
        i = tri[,1] - 1, j = tri[,2] - 1, k = tri[,3] - 1,
        type = "mesh3d",
        vertexcolor = hat_rgb,
        flatshading = TRUE,
        opacity = 1,
        name = paste0("phi_", i0),
        showlegend = FALSE
      )
      
      # ---- Hat wireframe edges ----
      x_hat_edges <- as.vector(t(cbind(x_hat[edges[,1]], x_hat[edges[,2]], NA)))
      y_hat_edges <- as.vector(t(cbind(y_hat[edges[,1]], y_hat[edges[,2]], NA)))
      z_hat_edges <- as.vector(t(cbind(z_hat[edges[,1]], z_hat[edges[,2]], NA)))
      
      p <- add_trace(p,
        x = x_hat_edges, y = y_hat_edges, z = z_hat_edges,
        type = "scatter3d", mode = "lines",
        line = list(color = hat_color[k], width = 2),
        showlegend = FALSE,
        hoverinfo = "none"
      )
    }
  }
  
   p <- add_trace(p,
    x = x_edges, y = y_edges, z = z_edges,
    type = "scatter3d", mode = "lines",
    line = list(color = "black", width = 3),
    showlegend = FALSE,
    hoverinfo = "none"
  )
   
  # ---- Nodes ----
  p <- add_trace(p,
    data = df3,
    x = ~x, y = ~y, z = ~z,
    mode = "markers",
    type = "scatter3d",
    marker = list(size = 4, color = "gray", symbol = 104),
    showlegend = FALSE
  )
  
  # ---- Scene ----
  p <- p %>% layout(font = list(family = "Palatino Linotype"),
    scene = global.scene.setter.rec.and.sphere(x_range, y_range, z_range)
  )
  
  return(p)
}


graph.plotter.3d.onecol <- function(graph, vec) {
  
  # Coordinates on the mesh
  x <- plotting.order(graph$mesh$V[, 1], graph)
  y <- plotting.order(graph$mesh$V[, 2], graph)
  z <- plotting.order(vec, graph)  # vector to plot
  
  # Axis ranges
  x_range <- range(x)
  y_range <- range(y)
  z_range <- range(z)
  z_range[1] <- z_range[1] - 10^-4
  
  # Vertical lines (vectorized)
  n <- length(x)
  x_vert <- rep(x, each = 3)
  y_vert <- rep(y, each = 3)
  z_vert <- as.vector(t(cbind(0, z, NA)))
  
  # Create plot
  p <- plot_ly() %>%
    
    # Main 3D curve over the graph
    add_trace(
      x = x, y = y, z = z,
      type = "scatter3d", mode = "lines",
      line = list(color = "blue", width = 3),
      showlegend = FALSE
    ) %>%
    # Graph base
    add_trace(
      x = x, y = y, z = z*0,
      type = "scatter3d", mode = "lines",
      line = list(color = "black", width = 3),
      showlegend = FALSE
    ) %>%
    
    # Vertical lines from base to curve
    add_trace(
      x = x_vert, y = y_vert, z = z_vert,
      type = "scatter3d", mode = "lines",
      line = list(color = "gray", width = 0.5),
      showlegend = FALSE
    ) %>%
    
    layout(
      scene = list(xaxis = list(title = "x", range = x_range),
              yaxis = list(title = "y", range = y_range),
              zaxis = list(
      title = "z",
      range = z_range,
      exponentformat = "e",
      tickformat = NULL#".4f"     # prevents SI prefixes like μ, m, k
    ),
              aspectratio = list(x = 2*(1+2/pi), 
                                 y = 2*(2/pi), 
                                 z = 1*(2/pi)),
              camera = list(eye = list(x = 5, 
                                       y = 3, 
                                       z = 3.5),
                            center = list(x = (1+2/pi)/2, 
                                          y = 0, 
                                          z = 0)))
    )
  
  return(p)
}


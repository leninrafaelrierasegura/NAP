## ------------------------------------------------------------------------------------------------------------------------------------------------
# library(ggplot2)
# library(reshape2)


## ------------------------------------------------------------------------------------------------------------------------------------------------
# loglog_line_equation <- function(x1, y1, slope) {
#   b <- log10(y1 / (x1 ^ slope))
# 
#   function(x) {
#     (x ^ slope) * (10 ^ b)
#   }
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# exp_line_equation <- function(x1, y1, slope) {
#   lnC <- log(y1) - slope * x1
# 
#   function(x) {
#     exp(lnC + slope * x)
#   }
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# compute_guiding_lines <- function(x_axis_vector, errors, theoretical_rates, line_equation_fun) {
#   guiding_lines <- matrix(NA, nrow = length(x_axis_vector), ncol = length(theoretical_rates))
# 
#   for (j in seq_along(theoretical_rates)) {
#     guiding_lines_aux <- matrix(NA, nrow = length(x_axis_vector), ncol = length(x_axis_vector))
# 
#     for (k in seq_along(x_axis_vector)) {
#       point_x1 <- x_axis_vector[k]
#       point_y1 <- errors[k, j]
#       slope <- theoretical_rates[j]
# 
#       line <- line_equation_fun(x1 = point_x1, y1 = point_y1, slope = slope)
#       guiding_lines_aux[, k] <- line(x_axis_vector)
#     }
# 
#     guiding_lines[, j] <- rowMeans(guiding_lines_aux)
#   }
# 
#   return(guiding_lines)
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# error.convergence.plotter <- function(x_axis_vector,
#                                       alpha_vector,
#                                       errors,
#                                       theoretical_rates,
#                                       observed_rates,
#                                       line_equation_fun,
#                                       fig_title,
#                                       x_axis_label,
#                                       apply_sqrt = FALSE) {
# 
#   x_vec <- if (apply_sqrt) sqrt(x_axis_vector) else x_axis_vector
# 
#   guiding_lines <- compute_guiding_lines(x_axis_vector = x_vec,
#                                          errors = errors,
#                                          theoretical_rates = theoretical_rates,
#                                          line_equation_fun = line_equation_fun)
# 
#   default_colors <- scales::hue_pal()(length(alpha_vector))
# 
#   plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
#     geom_line(
#       data = data.frame(x = x_vec, y = guiding_lines[, i]),
#       aes(x = x, y = y),
#       color = default_colors[i],
#       linetype = "dashed",
#       show.legend = FALSE
#     )
#   })
# 
#   df <- as.data.frame(cbind(x_vec, errors))
#   colnames(df) <- c("x_axis_vector", alpha_vector)
#   df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
# 
#   custom_labels <- paste0(formatC(alpha_vector, format = "f", digits = 2),
#                           " | ",
#                           formatC(theoretical_rates, format = "f", digits = 4),
#                           " | ",
#                           formatC(observed_rates, format = "f", digits = 4))
# 
#   df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)
# 
#   p <- ggplot() +
#     geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
#     geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
#     plot_lines +
#     labs(
#       title = fig_title,
#       x = x_axis_label,
#       y = expression(Error),
#       color = if (apply_sqrt) "alpha | bound benchmark | observed" else "alpha | theoretical | observed"
#     ) +
#     (if (apply_sqrt) {
#       scale_x_continuous(breaks = x_vec, labels = round(x_axis_vector, 4))
#     } else {
#       scale_x_log10(breaks = x_axis_vector, labels = round(x_axis_vector, 4))
#     }) +
#     (if (apply_sqrt) {
#       scale_y_continuous(trans = "log", labels = scales::scientific_format())
#     } else {
#       scale_y_log10(labels = scales::scientific_format())
#     }) +
#     theme_minimal() +
#     theme(text = element_text(family = "Palatino"),
#           legend.position = "bottom",
#           legend.direction = "vertical",
#           plot.margin = margin(0, 0, 0, 0),
#           plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
# 
#   return(p)
# }


## ------------------------------------------------------------------------------------------------------------------------------------------------
# errorConvergencePlotterLatex <- function(x_axis_vector,
#                                       alpha_vector,
#                                       errors,
#                                       theoretical_rates,
#                                       observed_rates,
#                                       line_equation_fun,
#                                       fig_title,
#                                       x_axis_label,
#                                       y_axis_label,
#                                       color_label,
#                                       apply_sqrt = FALSE) {
# 
#   relative_per_error <- 100 * abs(theoretical_rates - observed_rates) / abs(theoretical_rates)
# 
#   x_vec <- if (apply_sqrt) sqrt(x_axis_vector) else x_axis_vector
# 
#   guiding_lines <- compute_guiding_lines(x_axis_vector = x_vec,
#                                          errors = errors,
#                                          theoretical_rates = theoretical_rates,
#                                          line_equation_fun = line_equation_fun)
# 
#   default_colors <- scales::hue_pal()(length(alpha_vector))
# 
#   plot_lines <- lapply(1:ncol(guiding_lines), function(i) {
#     geom_line(
#       data = data.frame(x = x_vec, y = guiding_lines[, i]),
#       aes(x = x, y = y),
#       color = default_colors[i],
#       linetype = "dashed",
#       linewidth = 0.6,
#       show.legend = FALSE
#     )
#   })
# 
#   df <- as.data.frame(cbind(x_vec, errors))
#   colnames(df) <- c("x_axis_vector", alpha_vector)
#   df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
#   if (apply_sqrt) {
#     custom_labels <- sprintf("$%s \\; | \\; %s \\; | \\; %s$",
#       formatC(alpha_vector, format = "f", digits = 2),
#       formatC(theoretical_rates, format = "f", digits = 4),
#       formatC(observed_rates, format = "f", digits = 4))
#   } else {
#     custom_labels <- sprintf("$%s \\; | \\; %s \\; | \\; %s \\; | \\; %s$",
#       formatC(alpha_vector, format = "f", digits = 2),
#       formatC(theoretical_rates, format = "f", digits = 4),
#       formatC(observed_rates, format = "f", digits = 4),
#       formatC(relative_per_error, format = "f", digits = 2))
#   }
# 
#   df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)
# 
#   scientific_ticks <- function(x) {
#     labels <- formatC(x, format = "e", digits = 2)
#     labels <- sub("0+e", "e", labels)
#     labels <- sub("[.]e", "e", labels)
#     labels <- sub("e([+-])0+([0-9]+)$", "e\\1\\2", labels)
#     sub("e[+]", "e", labels)
#   }
# 
#   p <- ggplot() +
#     geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column), linewidth = 1.6) +
#     geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column), size = 2) +
#     plot_lines +
#     labs(
#       title = fig_title,
#       x = x_axis_label,
#       y = y_axis_label,
#       color = color_label
#     ) +
#     (if (apply_sqrt) {
#       scale_x_continuous(breaks = x_vec, labels = formatC(x_axis_vector, format = "f", digits = 0))
#     } else {
#       scale_x_log10(breaks = x_axis_vector, labels = scientific_ticks)
#     }) +
#     (if (apply_sqrt) {
#       #scale_y_continuous(trans = "log", labels = scales::scientific_format())
#       scale_y_continuous(trans = scales::log_trans(base = exp(1)), labels = scales::scientific_format())
#     } else {
#       scale_y_log10(labels = scales::scientific_format())
#     }) +
#     theme_minimal() +
#     # theme(text = element_text(family = "Palatino"),
#     #       legend.position = "bottom",
#     #       legend.direction = "vertical",
#     #       #plot.margin = margin(0, 0, 0, 0),
#     #       plot.title = element_text(hjust = 0.5, size = 18, face = "plain")) +
#     theme(
#   text = element_text(family = "Palatino", size = 18),
# 
#   plot.title = element_text(
#     hjust = 0.5,
#     size = 20,
#     face = "plain"
#   ),
# 
#   axis.title.x = element_text(size = 14),
#   axis.title.y = element_text(size = 14),
# 
#   axis.text.x = element_text(size = 12),
#   axis.text.y = element_text(size = 12),
# 
#   legend.title = element_text(size = 14),
#   legend.text = element_text(size = 12),
# 
#   legend.position = "bottom",
#   legend.direction = "vertical"
# )
# 
#   return(p)
# }


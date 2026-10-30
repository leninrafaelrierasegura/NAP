rough_convergence_plot <- function(result, corrected = FALSE) {
  p <- result$parameters
  errors <- if (corrected) result$corrected_errors else result$errors
  frame <- expand.grid(J = p$J, column = seq_along(p$alpha))
  frame$h <- 2^-frame$J
  frame$error <- as.vector(errors)
  frame$alpha <- factor(p$alpha[frame$column], levels = p$alpha)
  fit <- which(p$J %in% p$fit_J)
  guide <- frame
  for (a in seq_along(p$alpha)) {
    target <- guide$column == a
    model <- (2^-p$J)^p$alpha[a]
    if (!corrected) model <- model/p$J
    constant <- exp(mean(log(errors[fit, a]/model[fit])))
    guide$error[target] <- constant * model
  }
  label <- if (corrected) expression(J %.% E[h]) else expression(E[h])
  ggplot2::ggplot(frame, ggplot2::aes(h, error, color = alpha, group = alpha)) +
    ggplot2::geom_line(data = guide, linetype = "dashed", linewidth = .65, alpha = .8) +
    ggplot2::geom_line(linewidth = .8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_x_log10(breaks = 2^-p$J,
      labels = function(x) paste0("1/", round(1/x))) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_viridis_d(option = "D", begin = .05, end = .85) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2L, byrow = TRUE)) +
    ggplot2::labs(x = expression(h~"(log scale)"), y = label, color = expression(alpha),
      title = if (corrected) "Spatial convergence after logarithmic correction" else
        "Spatial convergence with rough forcing") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 13, face = "bold"))
}

rough_save_figures <- function(result, output_dir = file.path("data_files", "rough")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (corrected in c(FALSE, TRUE)) {
    plot <- rough_convergence_plot(result, corrected)
    name <- if (corrected) "conv_h_rough_corrected" else "conv_h_rough_raw"
    for (extension in c("png", "pdf"))
      ggplot2::ggsave(file.path(output_dir, paste0(name, ".", extension)),
        plot, width = 6.2, height = 4.7, dpi = 300, bg = "white")
  }
  invisible(output_dir)
}

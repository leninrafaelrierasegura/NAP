

error.convergence.plotter <- function(x_axis_vector, 
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
      show.legend = FALSE
    )
  })
  
  df <- as.data.frame(cbind(x_vec, errors))
  colnames(df) <- c("x_axis_vector", alpha_vector)
  df_melted <- melt(df, id.vars = "x_axis_vector", variable.name = "column", value.name = "value")
  
  custom_labels <- sprintf("$%s \\; | \\; %s \\; | \\; %s \\; | \\; %s$",
                           formatC(alpha_vector, format = "f", digits = 2),
                           formatC(theoretical_rates, format = "f", digits = 2),
                           formatC(observed_rates, format = "f", digits = 2),
                           formatC(relative_per_error, format = "f", digits = 2))
  
  df_melted$column <- factor(df_melted$column, levels = alpha_vector, labels = custom_labels)
  
  p <- ggplot() +
    geom_line(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    geom_point(data = df_melted, aes(x = x_axis_vector, y = value, color = column)) +
    plot_lines +
    labs(
      title = fig_title,
      x = x_axis_label,
      y = y_axis_label,
      color = color_label
    ) +
    (if (apply_sqrt) {
      scale_x_continuous(breaks = x_vec, labels = round(x_axis_vector, 4))
    } else {
      scale_x_log10(breaks = x_axis_vector, labels = round(x_axis_vector, 4))
    }) +
    (if (apply_sqrt) {
      #scale_y_continuous(trans = "log", labels = scales::scientific_format())
      scale_y_continuous(trans = scales::log_trans(base = exp(1)), labels = scales::scientific_format())
    } else {
      scale_y_log10(labels = scales::scientific_format())
    }) +
    theme_minimal() +
    theme(text = element_text(family = "Palatino"),
          legend.position = "right",
          legend.direction = "vertical",
          #plot.margin = margin(0, 0, 0, 0),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  
  return(p)
}




alpha_vector <- 0.4 + seq(0.1, 0.5, by = 0.05) + 2
m_vector <- c(1:14)
errors <- matrix(NA, nrow = length(m_vector), ncol = length(alpha_vector))
for (alpha in alpha_vector) {
  for (m in m_vector) {
    res <- rSPDE:::interp_rational_coefficients(
      order = m, 
      type_rational_approx = "brasil", 
      type_interp = "spline", 
      alpha = alpha)
    
    partFraction <- function(x) {
      term_matrix <- outer(x, res$p, "-") # matrix of (x - p[i])
      term_values <- sweep(1 / term_matrix, 2, res$r, "*") # multiply columns by r[i]
      return(rowSums(term_values) + res$k)
    }
    
    powerFun <- function(x){
      return(x^(-(alpha-floor(alpha))))
    }
    
    delta <- 10^(-(5+m)/2)
    x <- seq(1, 1/delta, length.out = 10)
    pf <- partFraction(x)
    po <- powerFun(x)
    Linferror <- max(abs(pf - po))
    
    m_index <- which(m_vector == m)
    alpha_index <- which(alpha_vector == alpha)
    errors[m_index, alpha_index] <- Linferror
  }
}

observed_rates <- rep(NA, length(alpha_vector))
for (u in 1:length(alpha_vector)) {
  observed_rates[u] <- coef(lm(log(errors[, u]) ~ sqrt(m_vector)))[2]
}

theoretical_rates <- - 2 * pi * sqrt(alpha_vector - floor(alpha_vector))

conv_rate_ratapprox <- error.convergence.plotter(x_axis_vector = m_vector, 
                                                 alpha_vector = alpha_vector - floor(alpha_vector), 
                                                 errors, 
                                                 theoretical_rates, 
                                                 observed_rates,
                                                 line_equation_fun = exp_line_equation,
                                                 fig_title = "",
                                                 x_axis_label = "$m\\quad(\\sqrt{\\phantom{x}}\\mbox{ scale})$",
                                                 y_axis_label = "$\\|f-r_m\\|_{L_{\\infty}([1,1/\\delta_m])}\\quad(\\log_{\\mathrm{e}}\\mbox{ scale})$",
                                                 color_label = "$\\quad\\;\\;\\{\\alpha\\}\\;|\\;\\mathtt{theo}\\;|\\;\\;\\mathtt{obs}\\;\\;|\\%\\mathtt{re}$",
                                                 apply_sqrt = TRUE)



# -----------------------------------------------------------------------------
# Convergence using either the generated CSV or rSPDE m1table,...,m4table.
# Computed coefficients use exact CSV beta values; package roots are interpolated.
# Run this section independently, or after the preceding package comparison.
# The CSV stores rho_m(t) = factor * prod(t-rc) / prod(t-rb), approximating
# t^(beta-1). Thus r_m(x) = rho_m(1/x)/x approximates x^(-beta), and equals
# factor * prod(1-rc*x) / prod(1-rb*x). Testing rho_m(x) against x^(-beta)
# directly would test the wrong rational function.
# -----------------------------------------------------------------------------
test_computed_fractional_roots <- function(
    csv_file = here::here("matlab", "computed_fractional_roots.csv"),
    grid_points = 10001L,
    coefficient_source = c("computed", "package"),
    interpolation = c("linear", "cubic"),
    beta_vector = NULL,
    max_m = NULL) {
  if (!is.null(max_m) &&
      (!is.numeric(max_m) || length(max_m) != 1L || !is.finite(max_m) ||
       max_m < 2 || max_m != floor(max_m) || max_m > .Machine$integer.max))
    stop("max_m must be NULL or an integer of at least 2 to fit a convergence rate.")
  coefficient_source <- match.arg(coefficient_source)
  interpolation <- match.arg(interpolation)
  if (coefficient_source == "package") {
    # Same tables as rSPDE:::m1table, ..., rSPDE:::m4table.
    if (is.null(beta_vector)) beta_vector <- c(.25, .45, .65)
    if (!is.numeric(beta_vector) || !length(beta_vector) ||
        any(!is.finite(beta_vector)) || any(beta_vector <= 0 | beta_vector >= 1) ||
        anyDuplicated(beta_vector)) stop("beta_vector must contain distinct values in (0, 1).")
    tab <- do.call(rbind, lapply(1:4, function(m) {
      mt <- get(paste0("m", m, "table"), envir = asNamespace("rSPDE"))
      if (any(beta_vector < min(mt$beta) | beta_vector > max(mt$beta)))
        stop("Requested beta is outside the package table range for m = ", m, ".")
      do.call(rbind, lapply(beta_vector, function(beta) {
        at <- function(column) {
          if (interpolation == "linear") approx(mt$beta, mt[[column]], xout = beta)$y else
            splinefun(mt$beta, mt[[column]], method = "fmm")(beta)
        }
        columns <- c("factor", if (m == 1L) "rc" else paste0("rc.", seq_len(m)),
                     paste0("rb.", seq_len(m+1L)))
        data.frame(m = m, beta = beta,
          kind = c("factor", rep("rc", m), rep("rb", m+1L)),
          index = c(0L, seq_len(m), seq_len(m+1L)),
          value = vapply(columns, at, numeric(1)), row.names = NULL)
      }))
    }))
  } else {
    tab <- read.csv(csv_file, stringsAsFactors = FALSE)
  }
  stopifnot(length(grid_points) == 1L, is.finite(grid_points), grid_points >= 3L,
            grid_points == as.integer(grid_points))
  numeric_columns <- c("m", "beta", "index", "value")
  if (!identical(names(tab), c("m", "beta", "kind", "index", "value")) ||
      !nrow(tab) || !all(vapply(tab[numeric_columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(tab[numeric_columns]))) ||
      any(tab$m < 1 | tab$m != round(tab$m)) ||
      any(tab$beta <= 0 | tab$beta >= 1) ||
      any(!tab$kind %in% c("factor", "rc", "rb")) ||
      anyDuplicated(tab[c("m", "beta", "kind", "index")]))
    stop("Malformed coefficient entries for source '", coefficient_source, "'.")
  if (coefficient_source == "computed" && !is.null(beta_vector)) {
    if (!is.numeric(beta_vector) || !length(beta_vector) ||
        any(!is.finite(beta_vector)) || anyDuplicated(beta_vector))
      stop("beta_vector must contain distinct finite values available in the CSV.")
    available <- sort(unique(tab$beta))
    selected <- vapply(beta_vector, function(beta) {
      index <- which(abs(available-beta) <= 32*.Machine$double.eps*max(1, abs(beta)))
      if (length(index) != 1L) stop("Computed coefficients are unavailable at beta = ", beta, ".")
      available[index]
    }, numeric(1))
    tab <- tab[tab$beta %in% selected, , drop = FALSE]
  }
  source_label <- if (coefficient_source == "package")
    paste0("rSPDE package coefficients (", interpolation, ")") else "Computed coefficients"

  available_orders <- sort(unique(tab$m))
  if (!is.null(max_m)) {
    if (max_m > max(available_orders))
      stop("Coefficient source '", coefficient_source, "' has a maximum available m of ",
           max(available_orders), "; requested max_m = ", max_m, ".")
    requested_orders <- seq_len(max_m)
    missing_orders <- setdiff(requested_orders, available_orders)
    if (length(missing_orders))
      stop("Missing coefficient orders: ", paste(missing_orders, collapse = ", "), ".")
    tab <- tab[tab$m %in% requested_orders, , drop = FALSE]
  }
  m_values <- sort(unique(tab$m))
  beta_values <- sort(unique(tab$beta))
  if (length(m_values) < 2L) stop("At least two orders are needed to assess convergence.")
  upper_limits <- 10^((5+m_values)/2)
  if (any(!is.finite(upper_limits))) stop("The requested sampling interval overflows.")
  common_upper <- min(upper_limits) # Shared subset of every generated interval.
  shape <- c(length(m_values), length(beta_values))
  labels <- list(m = m_values, beta = format(beta_values, trim = TRUE))
  moving_errors <- common_errors <- grid_changes <- matrix(
    NA_real_, shape[1], shape[2], dimnames = labels)

  for (j in seq_along(beta_values)) for (i in seq_along(m_values)) {
    m <- m_values[i]
    beta <- beta_values[j]
    rows <- tab[tab$m == m & tab$beta == beta, , drop = FALSE]
    values <- function(kind, indices) {
      part <- rows[rows$kind == kind, , drop = FALSE]
      if (nrow(part) != length(indices) || !setequal(part$index, indices))
        stop("Incomplete coefficients for m = ", m, ", beta = ", beta, ".")
      part$value[match(indices, part$index)]
    }
    factor <- values("factor", 0)
    rc <- values("rc", seq_len(m))
    rb <- values("rb", seq_len(m+1L))
    if (factor <= 0 || any(c(rc, rb) >= 0))
      stop("A positive factor and negative roots are required at m = ", m,
           ", beta = ", beta, ".")
    numerator <- sort(rc)
    denominator <- sort(rb)
    if (any(denominator[-length(denominator)] >= numerator) ||
        any(numerator >= denominator[-1L]))
      stop("Noninterlacing roots at m = ", m, ", beta = ", beta, ".")

    sampled_error <- function(upper, n) {
      # Logarithmic sampling resolves x close to 1 as well as the far tail.
      x <- exp(seq(0, log(upper), length.out = n))
      log_r <- log(factor)+rowSums(log1p(-outer(x, rc)))-
        rowSums(log1p(-outer(x, rb)))
      # r-x^(-beta) = x^(-beta)*expm1(log(r)+beta*log(x)).
      # This avoids subtracting nearly equal function values.
      max(abs(exp(-beta*log(x))*expm1(log_r+beta*log(x))))
    }
    for (domain in c("moving", "common")) {
      upper <- if (domain == "moving") upper_limits[i] else common_upper
      coarse <- sampled_error(upper, grid_points)
      refined <- sampled_error(upper, 2L*grid_points-1L)
      if (!is.finite(refined) || refined <= 0)
        stop("Invalid or unresolved approximation error at m = ", m, ", beta = ", beta, ".")
      if (domain == "moving") moving_errors[i, j] <- refined else
        common_errors[i, j] <- refined
      grid_changes[i, j] <- max(c(grid_changes[i, j],
        abs(refined-coarse)/refined), na.rm = TRUE)
    }
  }

  fit <- function(error) {
    model <- lm(log(error) ~ sqrt(m_values))
    c(slope = unname(coef(model)[2]), R_squared = summary(model)$r.squared)
  }
  moving_fits <- vapply(seq_along(beta_values), function(j) fit(moving_errors[, j]), numeric(2))
  common_fits <- vapply(seq_along(beta_values), function(j) fit(common_errors[, j]), numeric(2))
  # Same inverse-power benchmark as the preceding example. These Cheb-Pade
  # coefficients are not assumed to attain a best-approximation asymptotic rate.
  benchmark <- -2*pi*sqrt(beta_values)
  summary <- data.frame(beta = beta_values, benchmark = benchmark,
    moving_slope = moving_fits["slope", ], moving_R_squared = moving_fits["R_squared", ],
    moving_reduction = moving_errors[1, ]/moving_errors[nrow(moving_errors), ],
    moving_strictly_decreasing = vapply(seq_along(beta_values),
      function(j) all(diff(moving_errors[, j]) < 0), logical(1)),
    common_slope = common_fits["slope", ], common_R_squared = common_fits["R_squared", ],
    common_reduction = common_errors[1, ]/common_errors[nrow(common_errors), ],
    common_strictly_decreasing = vapply(seq_along(beta_values),
      function(j) all(diff(common_errors[, j]) < 0), logical(1)),
    max_sampling_relative_change = apply(grid_changes, 2, max), row.names = NULL)

  make_plot <- function(error, slopes, title, domain_label) {
    points <- do.call(rbind, lapply(seq_along(beta_values), function(j)
      data.frame(m = m_values, error = error[, j], beta = factor(beta_values[j], levels = beta_values))))
    guides <- do.call(rbind, lapply(seq_along(beta_values), function(j)
      data.frame(m = m_values, error = error[1, j]*
        exp(benchmark[j]*(sqrt(m_values)-sqrt(m_values[1]))),
        beta = factor(beta_values[j], levels = beta_values))))
    legend_labels <- sprintf("beta = %.2f | benchmark %.2f | fitted %.2f",
                              beta_values, benchmark, slopes)
    ggplot2::ggplot(points, ggplot2::aes(x = sqrt(m), y = error, colour = beta))+
      ggplot2::geom_line(linewidth = .8)+ggplot2::geom_point(size = 2)+
      ggplot2::geom_line(data = guides, linetype = "dashed", linewidth = .5)+
      ggplot2::scale_x_continuous(breaks = sqrt(m_values), labels = m_values)+
      ggplot2::scale_y_continuous(trans = scales::log_trans(base = exp(1)),
                                  labels = scales::label_scientific())+
      ggplot2::scale_colour_discrete(labels = legend_labels)+
      ggplot2::labs(title = title, subtitle = domain_label,
        x = "m (square-root scale)", y = "Sampled maximum absolute error (log scale)",
        colour = NULL, caption = "Dashed: inverse-power benchmark, not a required fitted slope.")+
      ggplot2::theme_minimal()+
      ggplot2::theme(legend.position = "bottom", legend.direction = "vertical",
                     plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  }
  if (any(summary$max_sampling_relative_change > .01))
    warning("Sampling refinement changed an error by more than 1%; increase grid_points.")
  if (any(!summary$moving_strictly_decreasing))
    message("Some errors are nonmonotone on the expanding intervals; inspect the common-interval comparison.")
  list(coefficient_source = coefficient_source, max_m = max(m_values),
    interpolation = if (coefficient_source == "package") interpolation else NULL,
    m_vector = m_values, beta_vector = beta_values, errors = moving_errors,
    errors_common = common_errors, summary = summary, grid_changes = grid_changes,
    common_interval = c(1, common_upper),
    plot = make_plot(moving_errors, moving_fits["slope", ],
      paste0(source_label, ": order-dependent intervals"), "x in [1, 1/delta_m]; delta_m = 10^(-(5+m)/2)"),
    plot_common = make_plot(common_errors, common_fits["slope", ],
      paste0(source_label, ": common interval"), paste0("x in [1, ", format(common_upper, scientific = TRUE), "] for every m")))
}

# Choose "computed" (CSV orders) or "package" (m=1:4).
coefficient_source <- "computed"
max_m <- 4L # Test m = 1,...,max_m; NULL uses all orders available in the source.
computed_coeff_test <- test_computed_fractional_roots(
  coefficient_source = coefficient_source,
  max_m = max_m,
  interpolation = "linear", # Package only: "linear" or "cubic".
  beta_vector = c(.25, .45, .65))
cat("\nCoefficient source:", computed_coeff_test$coefficient_source, "\n")
cat("\nSampled errors on [1, 1/delta_m]\n")
print(computed_coeff_test$errors)
cat("\nSampled errors on the common interval\n")
print(computed_coeff_test$errors_common)
cat("\nConvergence summary (benchmark comparison, not a correctness percentage)\n")
print(computed_coeff_test$summary)
conv_rate_computed_coefficients <- computed_coeff_test$plot
conv_rate_computed_coefficients_common <- computed_coeff_test$plot_common
# Keep plots in memory; no reports or experiment outputs are written here.
if (interactive()) {
  print(conv_rate_computed_coefficients)
  print(conv_rate_computed_coefficients_common)
}

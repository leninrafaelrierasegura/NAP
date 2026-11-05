alpha <- 0.35
m <- 2L

tt <- get(paste0("m_brasil", m, "t"), envir = asNamespace("rSPDE"))
brasil_row <- tt[abs(tt$alpha - alpha) < 1e-12, , drop = FALSE]
stopifnot(nrow(brasil_row) == 1L)

brasil <- list(
  r = as.numeric(unlist(brasil_row[paste0("r", seq_len(m))])),
  p = as.numeric(unlist(brasil_row[paste0("p", seq_len(m))])),
  k = brasil_row$k
)

wl2 <- rSPDE::rational.coefficients.wl2(
  alpha = alpha,
  d = 1,
  m = m,
  x_min = NULL,
  type = "covariance",
  k_term = TRUE
)

# Direct partial fraction: R(x) = k + sum_i r_i / (x - p_i).
partial_fraction <- function(x, r, p, k) {
  k + vapply(x, function(z) sum(r / (z - p)), numeric(1))
}

# With x = 1/lambda, R(1/x) = k + sum_i r_i*x / (1 - p_i*x).
# This evaluates the continuous extension at x = 0 without dividing by zero.
# Since alpha < 1, the WL2 integer factor is absent (q = 0).
reciprocal_partial_fraction <- function(x, r, p, k) {
  k + vapply(x, function(z) sum(r * z / (1 - p * z)), numeric(1))
}

x <- seq(0, 1, length.out = 1001L)
curves <- data.frame(
  x = x,
  brasil = partial_fraction(x, brasil$r, brasil$p, brasil$k),
  wl2 = partial_fraction(x, wl2$r, wl2$p, wl2$k),
  brasil_reciprocal = reciprocal_partial_fraction(x, brasil$r, brasil$p, brasil$k),
  wl2_reciprocal = reciprocal_partial_fraction(x, wl2$r, wl2$p, wl2$k),
  target = x^alpha
)

# Start the direct target where it equals the larger approximation at zero.
direct_target_max <- max(partial_fraction(0, brasil$r, brasil$p, brasil$k),
                         partial_fraction(0, wl2$r, wl2$p, wl2$k))
x_0 <- direct_target_max^(-1 / alpha)
stopifnot(is.finite(x_0), x_0 > 0, x_0 <= 1)
direct_target_x <- seq(x_0, 1, length.out = 1001L)

# Estimate maximum errors over the intervals where each target is plotted.
direct_error_x <- seq(x_0, 1, length.out = 100001L)
error_x <- seq(0, 1, length.out = 100001L)
max_errors <- data.frame(
  direct = vapply(list(brasil, wl2), function(cf) {
    max(abs(partial_fraction(direct_error_x, cf$r, cf$p, cf$k) -
              direct_error_x^(-alpha)))
  }, numeric(1)),
  reciprocal = vapply(list(brasil, wl2), function(cf) {
    max(abs(reciprocal_partial_fraction(error_x, cf$r, cf$p, cf$k) -
              error_x^alpha))
  }, numeric(1)),
  row.names = c("BRASIL", "Weighted L2")
)

plot_partial_fractions <- function() {
  old_par <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 5.2, 1))
  on.exit(par(old_par))
  colors <- c("#0072B2", "#D55E00")

  matplot(curves$x, curves[c("brasil", "wl2")], type = "l",
          col = colors, lty = c(1, 2), lwd = 2, xlim = c(0, 1),
          xlab = "x", ylab = "R(x)")
  lines(direct_target_x, direct_target_x^(-alpha),
        col = "#333333", lty = 3, lwd = 2)
  title(sprintf("Direct partial fractions (m = %d)", m),
        line = 3, cex.main = 1.1)
  mtext(expression(R(x) == k + sum(r[i] / (x - p[i]), i)),
        side = 3, line = 0.7, cex = 0.8)
  legend("topright",
         c(sprintf("%s: %.4g", rownames(max_errors), max_errors$direct),
           sprintf("Target: x^(-%g)", alpha)),
         title = sprintf("Grid max |error| on [x0, 1]\nx0 = %.4g", x_0),
         col = c(colors, "#333333"), lty = c(1, 2, 3), lwd = 2,
         bty = "n", cex = 0.85)

  matplot(curves$x, curves[c("brasil_reciprocal", "wl2_reciprocal", "target")],
          type = "l", col = c(colors, "#333333"), lty = c(1, 2, 3),
          lwd = 2, xlim = c(0, 1), xlab = "x", ylab = "R(1/x)")
  title(sprintf("Reciprocal-variable approximations (m = %d)", m),
        line = 3, cex.main = 1.1)
  mtext(expression(R(1/x) == k + sum(r[i] * x / (1 - p[i] * x), i)),
        side = 3, line = 0.7, cex = 0.8)
  legend("topleft",
         c(sprintf("%s: %.4g", rownames(max_errors), max_errors$reciprocal),
           sprintf("Target: x^%g", alpha)),
         title = "Grid max |error| (100,001 points)",
         col = c(colors, "#333333"), lty = c(1, 2, 3), lwd = 2,
         bty = "n", cex = 0.85)
}

plot_partial_fractions()


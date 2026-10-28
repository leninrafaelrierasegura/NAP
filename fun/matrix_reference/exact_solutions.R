## -----------------------------------------------------------------------------
# remotes::install_github("davidbolin/rspde", ref = "devel")
# remotes::install_github("davidbolin/metricgraph", ref = "devel")
library(rSPDE)
library(MetricGraph)
library(grateful)

library(ggplot2)
library(reshape2)
library(plotly)


## -----------------------------------------------------------------------------
# Functions to compute the exact solution to the fractional diffusion equation
g_linear <- function(r, A, lambda_j_alpha_half) {
  return(A * exp(-lambda_j_alpha_half * r))
  }
G_linear <- function(t, A) {
  return(A * t)
  }
g_exp <- function(r, A, mu) {
  return(A * exp(mu * r))
  }
G_exp <- function(t, A, lambda_j_alpha_half, mu) {
  exponent <- lambda_j_alpha_half + mu
  return(A * (exp(exponent * t) - 1) / exponent)
  }
g_poly <- function(r, A, n) {
  return(A * r^n)
}
G_poly <- function(t, A, lambda_j_alpha_half, n) {
  t <- as.vector(t)
  k_vals <- 0:n
  sum_term <- sapply(t, function(tt) {
    sum(((-lambda_j_alpha_half * tt)^k_vals) / factorial(k_vals))
  })
  coeff <- ((-1)^(n + 1)) * factorial(n) / (lambda_j_alpha_half^(n + 1))
  return(A * coeff * (1 - exp(lambda_j_alpha_half * t) * sum_term))
}
g_sin <- function(r, A, omega) {
  return(A * sin(omega * r))
}
G_sin <- function(t, A, lambda_j_alpha_half, omega) {
  denom <- lambda_j_alpha_half^2 + omega^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * sin(omega * t) - omega * cos(omega * t)) + omega
  return(A * numerator / denom)
}
g_cos <- function(r, A, theta) {
  return(A * cos(theta * r)) 
}
G_cos <- function(t, A, lambda_j_alpha_half, theta) {
  denom <- lambda_j_alpha_half^2 + theta^2
  numerator <- exp(lambda_j_alpha_half * t) * (lambda_j_alpha_half * cos(theta * t) + theta * sin(theta * t)) - lambda_j_alpha_half
  return(A * numerator / denom)
}


## -----------------------------------------------------------------------------
decayed_sine_convolution <- function(t, rate, A = 1, omega = pi) {
  stopifnot(all(is.finite(t)), all(t >= 0), length(rate) == 1L,
            is.finite(rate), rate >= 0, is.finite(A), is.finite(omega))
  z <- (rate + 1i * omega) * t
  small <- Mod(z) < 0.5
  answer <- numeric(length(t))
  if (any(small)) {
    # exp(-rate*t) * (exp(z)-1)/z; subtract the real constant before
    # taking the imaginary part so the sine response is O(t^2) accurately.
    zz <- z[small]
    term <- rep(1 + 0i, length(zz))
    series <- rep(0 + 0i, length(zz))
    for (k in 1:24) {
      term <- term * zz / (k + 1)
      series <- series + term
    }
    answer[small] <- A * exp(-rate * t[small]) * t[small] * Im(series)
  }
  if (any(!small)) {
    tt <- t[!small]
    answer[!small] <- A * Im(
      (exp(1i * omega * tt) - exp(-rate * tt)) / (rate + 1i * omega)
    )
  }
  answer
}

exact_modal_coefficients <- function(t, rates, initial, forcing, A = 1, omega = pi) {
  stopifnot(length(rates) == length(initial), length(rates) == length(forcing),
            all(is.finite(rates)), all(rates > 0))
  out <- matrix(0, nrow = length(rates), ncol = length(t))
  for (j in seq_along(rates)) {
    out[j, ] <- initial[j] * exp(-rates[j] * t) + forcing[j] *
      decayed_sine_convolution(t, rates[j], A = A, omega = omega)
  }
  out
}

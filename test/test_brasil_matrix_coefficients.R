# Run from the project root in its R environment. No saved results or messages
# to external services are produced.
library(Matrix)
source("fun/brasil_coefficients.R")

orders <- brasil_supported_orders()
betas <- c(seq(0.25, 0.95, 0.1), 1/3)
worst_reconstruction <- 0
valid_cases <- character()
rejected_cases <- 0L
for (m in orders) for (beta in betas) {
  # Compare extraction against rSPDE's source columns, independently of H.
  table <- get(paste0("m_brasil", m, "t"), asNamespace("rSPDE"))
  column <- function(name) stats::approx(table$alpha, table[[name]], beta)$y
  raw_poles <- vapply(paste0("p", seq_len(m)), column, numeric(1))
  raw_residues <- vapply(paste0("r", seq_len(m)), column, numeric(1))
  # Some installed high-order tables do not have the positive structure.
  # Verify explicit rejection instead of clipping or altering their values.
  if (any(raw_poles >= 0) || any(raw_residues <= 0) || column("k") < 0 ||
      anyDuplicated(raw_poles)) {
    stopifnot(inherits(try(brasil_coefficients(m, beta), silent = TRUE), "try-error"))
    rejected_cases <- rejected_cases + 1L
    next
  }
  coefficients <- brasil_coefficients(m, beta)
  valid_cases <- c(valid_cases, paste(m, beta))
  y <- 10^seq(0, 12, length.out = 401L)
  expected <- rep(column("k"), length(y))
  for (j in seq_len(m)) expected <- expected + column(paste0("r", j)) /
    (y - column(paste0("p", j)))
  stopifnot(max(abs(brasil_inverse(y, coefficients) - expected)) < 1e-12,
            max(abs(brasil_power(1/y, coefficients) - expected)) < 1e-12,
            identical(brasil_power(0, coefficients), coefficients$constant))
  for (strength in c(2e-8, 1e-4, 0.1, 10)) {
    transfer <- brasil_partial_fractions(coefficients, strength, scaling = 1)
    stopifnot(length(transfer$p) == m, !is.complex(transfer$p),
              all(transfer$p < coefficients$poles),
              all(transfer$p[-1L] > head(coefficients$poles, -1L)),
              all(transfer$r > 0), transfer$k >= 0, transfer$k < 1)
    actual <- transfer$k + as.numeric((1 / outer(y, transfer$p, "-")) %*% transfer$r)
    gap <- max(abs(actual - expected/(expected + strength)))
    worst_reconstruction <- max(worst_reconstruction, gap)
    stopifnot(gap < 2e-12)
  }
}
cat("BRASIL extraction, reciprocal evaluation and transfer checks:",
    length(valid_cases) * 4L, "cases; max discrepancy",
    format(worst_reconstruction, scientific = TRUE), "\n")
cat("Invalid signed table rows explicitly rejected:", rejected_cases, "\n")

isolated <- new.env()
source("fun/isolated_rational_convergence.R", local = isolated)
matrix_reference <- new.env()
source("fun/matrix_reference/fractional_evolution.R", local = matrix_reference)
worst_matrix_gap <- 0
for (h in c(0.25, 0.1)) {
  data <- isolated$isolated_tadpole_data(h, kappa = 4)
  spectrum <- isolated$isolated_discrete_spectrum(data$C, data$L)
  rhs <- as.numeric(data$C %*% data$initial + 0.03 * data$spatial_load)
  for (helper in list(matrix_reference, isolated)) {
    for (beta in seq(0.25, 0.95, 0.1)) for (m in intersect(c(1L, 4L, 14L, 20L), orders)) {
      if (!paste(m, beta) %in% valid_cases) next
      for (tau in c(2e-5, 0.2)) for (diffusion in c(1e-3, 1)) {
        op <- helper$my.fractional.operators.frac(data$L, beta, data$C,
          scale.factor = data$kappa^2, m = m, time_step = tau, diffusion = diffusion)
        coefficients <- brasil_coefficients(m, beta)
        inverse <- brasil_inverse(spectrum$values/data$kappa^2, coefficients)
        multiplier <- inverse/(inverse + tau * diffusion * data$kappa^(2*beta))
        expected <- as.numeric(spectrum$vectors %*%
          (multiplier * as.numeric(crossprod(spectrum$vectors, rhs))))
        actual <- as.numeric(helper$my.solver.frac(op, rhs))
        gap <- max(abs(actual - expected))
        worst_matrix_gap <- max(worst_matrix_gap, gap)
        stopifnot(length(op$residues) == m, length(op$factorizations) == m,
                  op$remainder > 0, !is.null(op$mass_factor),
                  gap < 2e-10 * max(1, max(abs(expected))))
      }
    }
    op <- helper$my.fractional.operators.frac(data$L, 1, data$C,
      scale.factor = data$kappa^2, m = 4L, time_step = 0.03, diffusion = 0.7)
    expected <- as.numeric(solve(data$C + 0.03 * 0.7 * data$L, rhs))
    stopifnot(max(abs(as.numeric(helper$my.solver.frac(op, rhs)) - expected)) < 1e-11)
  }
}
cat("Both sparse matrix helpers agree with independent complete-eigensystem",
    "evaluation; max discrepancy", format(worst_matrix_gap, scientific = TRUE), "\n")

# Exercise complete matrix time histories and the isolated study's existing
# independent scalar diagnostic, including a non-tabulated exponent.
time_seq <- seq(0, 1, length.out = 51L)
results <- list()
for (alpha in c(seq(0.5, 1.9, 0.2), 2/3, 2)) {
  reference <- isolated$isolated_exact_be_reference(data$C, spectrum,
    data$initial, data$spatial_load, time_seq, alpha, diffusion = 0.3)
  for (m in 1:4) {
    result <- isolated$isolated_rational_error(data, spectrum, reference, m)
    stopifnot(is.finite(result$error), result$absolute_error_disagreement < 1e-9,
              result$final_state_gap_C < 1e-9,
              identical(result$coefficient_source, "brasil"))
    if (alpha == 2) stopifnot(result$error < 1e-10)
    results[[length(results) + 1L]] <- result
  }
}
cat("Isolated time-history checks:", length(results), "cases passed.\n")

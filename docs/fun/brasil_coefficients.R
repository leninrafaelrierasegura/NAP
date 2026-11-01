# Coefficients only: the studies retain their assembled FEM matrix solvers.
brasil_supported_orders <- function() {
  tables <- grep("^m_brasil[0-9]+t$", ls(asNamespace("rSPDE")), value = TRUE)
  sort(as.integer(sub("^m_brasil([0-9]+)t$", "\\1", tables)))
}

brasil_coefficients <- function(m, beta) {
  stopifnot(length(m) == 1L, is.finite(m), m %in% brasil_supported_orders(),
            length(beta) == 1L, is.finite(beta), beta > 0, beta < 1)
  table <- get(paste0("m_brasil", m, "t"), envir = asNamespace("rSPDE"))
  tolerance <- 32 * .Machine$double.eps * max(1, abs(beta))
  index <- which(abs(table$alpha - beta) <= tolerance)
  value <- function(column) {
    if (length(index) == 1L) table[[column]][index] else
      stats::approx(table$alpha, table[[column]], xout = beta)$y
  }
  poles <- vapply(paste0("p", seq_len(m)), value, numeric(1))
  residues <- vapply(paste0("r", seq_len(m)), value, numeric(1))
  constant <- value("k")
  if (any(!is.finite(c(poles, residues, constant))) || any(poles >= 0) ||
      any(residues <= 0) || constant < 0 || anyDuplicated(poles))
    stop("BRASIL table at m = ", m, ", beta = ", beta,
         " requires distinct negative poles, positive residues, and a nonnegative constant.")
  index <- order(poles)
  list(poles = unname(poles[index]), residues = unname(residues[index]),
       constant = unname(constant), m = as.integer(m), beta = beta)
}

# S(y) approximates y^(-beta); B(x) = S(1/x) approximates x^beta.
brasil_inverse <- function(y, coefficients) {
  stopifnot(all(is.finite(y)), all(y > 0))
  coefficients$constant + as.numeric(
    (1 / outer(y, coefficients$poles, "-")) %*% coefficients$residues)
}

brasil_power <- function(x, coefficients) {
  stopifnot(all(is.finite(x)), all(x >= 0))
  coefficients$constant + as.numeric(
    (x / (1 - outer(x, coefficients$poles))) %*% coefficients$residues)
}

# Partial fractions of H(y) = S(y)/(S(y) + time_step * scaling).
# Bracket the combined denominator's real roots between the S poles,
# avoiding ill-conditioned polynomial expansion and complex computed roots.
brasil_partial_fractions <- function(coefficients, time_step, scaling) {
  strength <- time_step * scaling
  stopifnot(length(strength) == 1L, is.finite(strength), strength > 0)
  p <- coefficients$poles
  r <- coefficients$residues
  constant <- coefficients$constant
  poles <- numeric(length(p))
  denominator <- function(y) constant + strength + sum(r / (y - p))
  for (i in seq_along(p)) {
    upper <- p[i]
    lower <- if (i == 1L) 2 * upper else p[i - 1L]
    if (i == 1L) while (denominator(lower) <= 0) lower <- 2 * lower
    width <- upper - lower
    normalized <- function(t) {
      if (t == 0 && i > 1L) return(.Machine$double.xmax)
      if (t == 1) return(-.Machine$double.xmax)
      denominator(lower + width * t)
    }
    root <- stats::uniroot(normalized, c(0, 1),
      tol = 2 * .Machine$double.eps, maxiter = 1000L)$root
    poles[i] <- lower + width * root
  }
  residues <- vapply(poles, function(z) strength / sum(r / (z - p)^2), numeric(1))
  remainder <- constant / (constant + strength)
  stopifnot(all(is.finite(c(poles, residues, remainder))), all(poles < 0),
            all(residues > 0), remainder >= 0, remainder < 1)
  y <- 10^seq(0, 12, length.out = 257L)
  inverse <- brasil_inverse(y, coefficients)
  expected <- inverse / (inverse + strength)
  actual <- remainder + as.numeric((1 / outer(y, poles, "-")) %*% residues)
  reconstruction_error <- max(abs(actual - expected))
  if (reconstruction_error > 2e-12)
    stop("Inaccurate BRASIL time-step partial fractions: ", reconstruction_error)
  list(p = poles, r = residues, k = remainder,
       reconstruction_error = reconstruction_error)
}
